using Pkg
Pkg.activate(".")
Pkg.instantiate()

DEPOT_PATH

using NumericalEarth, Oceananigans, Oceananigans.Units
using Oceananigans.BoundaryConditions
using Oceananigans.ImmersedBoundaries: InterfaceImmersedCondition
using Oceananigans.Fields: interpolate!
using Oceananigans.Operators: Δzᶠᶜᶜ, Δzᶜᶠᶜ
using Oceananigans.ImmersedBoundaries: immersed_peripheral_node, immersed_inactive_node
using Oceananigans.OrthogonalSphericalShellGrids: RotatedLatitudeLongitudeGrid
using Oceananigans.TurbulenceClosures: VerticalFormulation
using NumericalEarth.DataWrangling: SurfaceFluxRestoring
using Dates, Printf, Base64
using Oceananigans.OutputReaders: Clamp, Linear, Cyclical
using CopernicusMarine  

# Oceananigans.defaults.FloatType = Float32

using CUDA
arch = GPU()

# using Metal
# Oceananigans.defaults.FloatType = Float32
# metal = Metal.MetalBackend()
# arch = GPU(metal)

# arch = CPU() ########### Change to GPU() to run on a GPU

mp4_html(path) = HTML(string(""))

# =====================================================================
#  REGION CONFIG.  Edit ONLY this cell to move the model.
# =====================================================================

# Are you focused on a POLE / high latitudes? 
# true  → rotated grid (moves the geographic pole out of the domain so
#         cells don't shrink to a singularity). Needed above ~75° or if
#         the domain contains a pole.
# false → plain lat-lon grid (correct and simpler for mid/low latitude).
const POLAR = false
 
# Region bounds (real geographic degrees; lon −180..180)
const λ₁, λ₂ = -50, 10     # west, east longitude
const φ₁, φ₂ = -10.0, 10.0     # south, north latitude

# Coordinates for the Arctic, with POLAR = true
#const λ₁, λ₂ = 151, 212   # rotated longitude (east ↔ west edges → continents)
#const φ₁, φ₂ = -34.8, 27  # rotated latitude: south = Atlantic (~50°N) ↔ north (cuts across Bering Sea)
 
# Rotated-pole location (ONLY used when POLAR = true)
# Choose a north_pole (λ, φ) that places the real pole OUTSIDE your
# domain. For the Arctic cap the spinup used (145, 0).
const NORTH_POLE = (145, 0)
 
# Resolution & vertical grid
const RES   = 1/12        # degrees per cell 1 = 1degree resolution, 1/4 = quarter degree resolution etc.
const DEPTH = 5000meters      # ≥ deepest point in your region
  
# Forcing time window
const DATE_START = DateTime(2010, 9, 1)
const DATE_STOP  = DateTime(2010, 11, 1)
 
# A human label for plot titles / filenames
const REGION_NAME = "Equator" # This is used in the files that are saved, so change to the region you are looking at
 
@info "region $REGION_NAME: lon ($λ₁, $λ₂) lat ($φ₁, $φ₂), res $(RES)°, polar=$POLAR" 

# =====================================================================
#  GRID.  Rotated or plain, chosen by POLAR. No edits needed.
# =====================================================================

using Oceananigans.OrthogonalSphericalShellGrids: RotatedLatitudeLongitudeGrid
 
Nx = round(Int, (λ₂ - λ₁) / RES)
Ny = round(Int, (φ₂ - φ₁) / RES)

# vertical grid: constant spacing at the surface and bottom, stretched coordinates in between.
z = ReferenceToStretchedDiscretization(; extent = DEPTH, bias=:right, bias_edge=0, 
    constant_spacing = 20, constant_spacing_extent = 300, 
    stretching = PowerLawStretching(1.02), maximum_spacing = 150, mutable = false)

Nz = length(z) # number of vertical levels

if POLAR
    underlying_grid = RotatedLatitudeLongitudeGrid(arch;
        size       = (Nx, Ny, Nz),
        longitude  = (λ₁, λ₂),
        latitude   = (φ₁, φ₂),
        north_pole = NORTH_POLE,
        z,
        halo = (7, 7, 7))
else
    underlying_grid = LatitudeLongitudeGrid(arch;
        size      = (Nx, Ny, Nz),
        longitude = (λ₁, λ₂),
        latitude  = (φ₁, φ₂),
        z,
        topology = (Bounded, Bounded, Bounded), 
        halo = (7, 7, 7))

    GEBCO_grid = LatitudeLongitudeGrid(arch;
        size      = (Nx, Ny),
        longitude = (λ₁, λ₂),
        latitude  = (φ₁, φ₂),
        halo = (7, 7), 
        topology = (Bounded, Bounded, Flat))
end

@info "grid $(Nx)×$(Ny)×$(Nz)  ($(POLAR ? "rotated" : "plain") lat-lon)"

# fill in the below with your credentials for NASA Earth Data and Copernicus Marine Data Store
ENV["ECCO_USERNAME"] = "nicolaos.petropoulos"
ENV["ECCO_WEBDAV_PASSWORD"] = "ED-Azer-ty;09!"
ENV["COPERNICUS_USERNAME"] = "npetropoulos96@gmail.com"
ENV["COPERNICUS_PASSWORD"] = "Copernicus-Azer-ty;09!"

# Set these to the directories you want the forcing data to download to: GLORYS is the boundary forcing
# and JRA55 is the atmospheric forcing
ENV["GLORYS_DIR"] = "Forcings/GLORYSDaily"
ENV["JRA55_DIR"] = "Forcings/JRA55"

dir_kw       = haskey(ENV, "DATA_DIR")   ? (; dir = ENV["DATA_DIR"])   : (;)
dir_kwGLORYS = haskey(ENV, "GLORYS_DIR") ? (; dir = ENV["GLORYS_DIR"]) : (;)
dir_kwJRA55  = haskey(ENV, "JRA55_DIR")  ? (; dir = ENV["JRA55_DIR"])  : (;)

# =====================================================================
#  BATHYMETRY
# =====================================================================

# ETOPO2022 bathymetry 
#= bathymetry = Metadatum(:bottom_height; dataset = ETOPO2022(), dir_kw...)
bottom_height = regrid_bathymetry(underlying_grid, bathymetry;
                                  minimum_depth = 15,
                                  interpolation_passes = 2,
                                  major_basins = 1)   # keep the single connected basin
 
grid = ImmersedBoundaryGrid(underlying_grid, PartialCellBottom(bottom_height); active_cells_map = true) =#

# Custom bathymetry. Here, we use GEBCO2026, downloaded from https://download.gebco.net/
# and saved in Topography/ as a .nc file
λ = λnodes(underlying_grid, Center())
φ = φnodes(underlying_grid, Center())

using NCDatasets, Interpolations

ds = Dataset("Topography/topography.nc") # Replace with your downloaded file name
raw_lon = ds["lon"][:]
raw_lat = ds["lat"][:]
raw_elevation = ds["elevation"][:, :]

gebco_interp = LinearInterpolation((raw_lon, raw_lat), raw_elevation)
bathymetry_matrix = [gebco_interp(l, p) for l in λ, p in φ]

elevation_field = CenterField(GEBCO_grid)
set!(elevation_field, bathymetry_matrix)

bottom_height = Field{Center, Center, Nothing}(underlying_grid)
interpolate!(bottom_height, elevation_field)

grid = ImmersedBoundaryGrid(underlying_grid, PartialCellBottom(bottom_height); active_cells_map = true)
@info grid 

# =====================================================================
#  PLOT the bathymetry / land mask to check the region.
# =====================================================================
using CairoMakie, SixelTerm
 
h_bottom = Array(interior(bottom_height_field(grid), :, :, 1))
h_bottom[h_bottom .≥ 0] .= NaN
 
fig = Figure(size = (700, 600))
ax = Axis(fig[1, 1],
          xlabel = POLAR ? "grid i (rotated)" : "longitude [°]",
          ylabel = POLAR ? "grid j (rotated)" : "latitude [°]",
          title  = "$REGION_NAME bathymetry ($(POLAR ? "rotated" : "plain") grid)",
          aspect = DataAspect())
hm = heatmap!(ax, λ, φ, h_bottom, colormap = :deep, colorrange = (-DEPTH, 0), nan_color = :gray80)
Colorbar(fig[1, 2], hm, label = "bottom height [m]")
save("bathymetry.pdf", fig, pdf_version="1.4")

# Which edges are OPEN ocean boundaries?
# The open boundaries will be forced with the GLORYSDaily dataset (daily temperature, salinity and velocities)
# true  = open (needs boundary forcing); false = closed coastline.
const OPEN = (west = true, east = true, south = true, north = true)

# =====================================================================
#  BOUNDARY CONDITIONS — built automatically from the OPEN config.
#
#  For every edge marked `true` in OPEN, this builds:
#    - a tracer OBC (T, S)  fed by GLORYS at that edge
#    - a normal-velocity OBC (u on west/east, v on south/north)
#    - a barotropic Flather OBC on the normal transport
#  and a sponge mask that relaxes inward from ALL open edges.
#
#  Edges NOT in OPEN (or false) are left closed (solid coastline).
#
#  Indexing / staggering per edge:
#     west  : i = 1     normal = u   Δzᶠᶜᶜ  loc (Face,   Center)
#     east  : i = Nx    normal = u   Δzᶠᶜᶜ  loc (Face,   Center)
#     south : j = 1     normal = v   Δzᶜᶠᶜ  loc (Center, Face)
#     north : j = Ny    normal = v   Δzᶜᶠᶜ  loc (Center, Face)
# =====================================================================

# forcing time series (GLORYS), region tracks the domain
dates   = DATE_START : Day(1) : DATE_STOP
dataset = GLORYSDaily()
# pad the forcing region slightly beyond the model bounds so edges are covered
region  = BoundingBox(longitude = (λ₁ - 1, λ₂ + 1), latitude = (φ₁ - 1, φ₂ + 1))

Tᵉˣᵗ = FieldTimeSeries(Metadata(:temperature; dates, dataset, region, dir_kwGLORYS...), grid, inpainting=100, time_indices_in_memory=2, time_indexing=Clamp())
Sᵉˣᵗ = FieldTimeSeries(Metadata(:salinity; dates, dataset, region, dir_kwGLORYS...), grid, inpainting=100, time_indices_in_memory=2, time_indexing=Clamp())
uᵉˣᵗ = FieldTimeSeries(Metadata(:u_velocity; dates, dataset, region, dir_kwGLORYS...), grid, inpainting=100, time_indices_in_memory=2, time_indexing=Clamp())
vᵉˣᵗ = FieldTimeSeries(Metadata(:v_velocity; dates, dataset, region, dir_kwGLORYS...), grid, inpainting=100, time_indices_in_memory=2, time_indexing=Clamp())
ηᵉˣᵗ = FieldTimeSeries(Metadata(:free_surface; dates, dataset, region, dir_kwGLORYS...), grid, inpainting=100, time_indices_in_memory=2, time_indexing=Clamp())
nothing #hide

@inline wetcell(i, j, k, grid, ℓx, ℓy, ℓz) =
    !immersed_peripheral_node(i, j, k, grid, ℓx, ℓy, ℓz) & !immersed_inactive_node(i, j, k, grid, ℓx, ℓy, ℓz)

@inline function vertical_integral(i, j, grid, vel, t, Δz, ℓx, ℓy, ℓz)
    U = zero(eltype(grid))
    @inbounds for k in 1:grid.Nz
        wet = wetcell(i, j, k, grid, ℓx, ℓy, ℓz)
        U += ifelse(wet, vel[i, j, k, t] * Δz(i, j, k, grid), zero(U))
    end
    return U
end

# tracer boundary functions (one per edge)
# west/east are i-edges → the free index is j; south/north are j-edges → free index is i.
@inline west_trc(j, k, grid, clock, fields, φ) = @inbounds φ[1, j, k, Oceananigans.Units.Time(clock.time)]
@inline east_trc(j, k, grid, clock, fields, φ) = @inbounds φ[grid.Nx, j, k, Oceananigans.Units.Time(clock.time)]
@inline south_trc(i, k, grid, clock, fields, φ) = @inbounds φ[i, 1, k, Oceananigans.Units.Time(clock.time)]
@inline north_trc(i, k, grid, clock, fields, φ) = @inbounds φ[i, grid.Ny, k, Oceananigans.Units.Time(clock.time)]

# normal-velocity boundary functions 
# west/east: normal component is u, indexed at the i-face
@inline west_u(j, k, grid, clock, fields, φ) = @inbounds φ[1, j, k, Oceananigans.Units.Time(clock.time)]
@inline east_u(j, k, grid, clock, fields, φ) = @inbounds φ[grid.Nx+1, j, k, Oceananigans.Units.Time(clock.time)]
# south/north: normal component is v, indexed at the j-face
@inline south_v(i, k, grid, clock, fields, φ) = @inbounds φ[i, 1, k, Oceananigans.Units.Time(clock.time)]
@inline north_v(i, k, grid, clock, fields, φ) = @inbounds φ[i, grid.Ny+1, k, Oceananigans.Units.Time(clock.time)]

# assemble OBC keyword pairs only for the OPEN edges
rad = NormalRadiation(inflow_timescale = 1days)
# rad = PerturbationAdvection(inflow_timescale = 1days)
# rad = nothing

T_pairs = Pair{Symbol,Any}[]; S_pairs = Pair{Symbol,Any}[]
u_pairs = Pair{Symbol,Any}[]; v_pairs = Pair{Symbol,Any}[]


if get(OPEN, :west, false)
    push!(T_pairs, :west => ValueBoundaryCondition(west_trc, discrete_form=true, parameters=Tᵉˣᵗ, scheme=rad))
    push!(S_pairs, :west => ValueBoundaryCondition(west_trc, discrete_form=true, parameters=Sᵉˣᵗ, scheme=rad))
    push!(u_pairs, :west => NormalFlowBoundaryCondition(west_u, discrete_form=true, parameters=uᵉˣᵗ, scheme=rad))
end
if get(OPEN, :east, false)
    push!(T_pairs, :east => ValueBoundaryCondition(east_trc, discrete_form=true, parameters=Tᵉˣᵗ, scheme=rad))
    push!(S_pairs, :east => ValueBoundaryCondition(east_trc, discrete_form=true, parameters=Sᵉˣᵗ, scheme=rad))
    push!(u_pairs, :east => NormalFlowBoundaryCondition(east_u, discrete_form=true, parameters=uᵉˣᵗ, scheme=rad))
end
if get(OPEN, :south, false)
    push!(T_pairs, :south => ValueBoundaryCondition(south_trc, discrete_form=true, parameters=Tᵉˣᵗ, scheme=rad))
    push!(S_pairs, :south => ValueBoundaryCondition(south_trc, discrete_form=true, parameters=Sᵉˣᵗ, scheme=rad))
    push!(v_pairs, :south => NormalFlowBoundaryCondition(south_v, discrete_form=true, parameters=vᵉˣᵗ, scheme=rad))
end
if get(OPEN, :north, false)
    push!(T_pairs, :north => ValueBoundaryCondition(north_trc, discrete_form=true, parameters=Tᵉˣᵗ, scheme=rad))
    push!(S_pairs, :north => ValueBoundaryCondition(north_trc, discrete_form=true, parameters=Sᵉˣᵗ, scheme=rad))
    push!(v_pairs, :north => NormalFlowBoundaryCondition(north_v, discrete_form=true, parameters=vᵉˣᵗ, scheme=rad))
end

T_obcs = FieldBoundaryConditions(; T_pairs...)
S_obcs = FieldBoundaryConditions(; S_pairs...)
u_obcs = FieldBoundaryConditions(; u_pairs...)
v_obcs = FieldBoundaryConditions(; v_pairs...)

# sponge layer at open boundaries
sponge_width = 4   # cells
sponge_mask  = zeros(Nx, Ny, Nz)
for j in 1:Ny, i in 1:Nx
    dists = Float64[]
    get(OPEN, :west,  false) && push!(dists, i - 1)
    get(OPEN, :east,  false) && push!(dists, Nx - i)
    get(OPEN, :south, false) && push!(dists, j - 1)
    get(OPEN, :north, false) && push!(dists, Ny - j)
    d = isempty(dists) ? Inf : minimum(dists)      # nearest open edge
    sponge_mask[i, j, :] .= exp(-(d / sponge_width)^2)
end

sponge_width_TS = 4   # cells
sponge_mask_TS  = zeros(Nx, Ny, Nz)
for j in 1:Ny, i in 1:Nx, k in 1:Nz
    dists = Float64[]
    get(OPEN, :west,  false) && push!(dists, i - 1)
    get(OPEN, :east,  false) && push!(dists, Nx - i)
    get(OPEN, :south, false) && push!(dists, j - 1)
    get(OPEN, :north, false) && push!(dists, Ny - j)
    # push!(dists, Nz - k)
    d = isempty(dists) ? Inf : minimum(dists)      # nearest open edge
    if d < Nz - k
        sponge_mask_TS[i, j, k] = exp(-(d / sponge_width)^2)
    else
        sponge_mask_TS[i, j, k] = exp(-((Nz - k) / sponge_width_TS)^2)
    end
end

FT = DatasetRestoring(Metadata(:temperature; dates, dataset, region, dir_kwGLORYS...), grid; rate = 1/1days, inpainting=100, time_indexing=Clamp())
Fu = DatasetRestoring(Metadata(:u_velocity;  dates, dataset, region, dir_kwGLORYS...), grid; rate = 1/20minutes, inpainting=100, time_indexing=Clamp())
Fv = DatasetRestoring(Metadata(:v_velocity;  dates, dataset, region, dir_kwGLORYS...), grid; rate = 1/20minutes, inpainting=100, time_indexing=Clamp())
FS = DatasetRestoring(Metadata(:salinity;    dates, dataset, region, dir_kwGLORYS...), grid; rate = 1/1days, inpainting=100, time_indexing=Clamp())

# Tides
harmonics = TidalHarmonics(constituents = (:M2,),
                           frequencies = (1.405189e-4,),
                           phases = (1.7,),
                           equilibrium_amplitudes = (0.168,),
                           species = (2,))

FTides = tidal_forcing(harmonics)
# BCTides = tidal_boundary_conditions(grid, harmonics)

@info "open edges: $(join([string(e) for (e,v) in pairs(OPEN) if v], ", "))"

# Turbulent closure 
closure = CATKEVerticalDiffusivity(minimum_tke=1e-7)
# closure = ScalarDiffusivity(ExplicitTimeDiscretization(), VerticalFormulation(), Float64; ν = 1e-5, κ = 1e-5)
time_discretization = AdaptiveVerticallyImplicitDiscretization(cfl=0.5)

# =====================================================================
#  OCEAN MODEL
# =====================================================================

ocean = ocean_simulation(grid;
                         free_surface = SplitExplicitFreeSurface(grid; substeps=30),
                         momentum_advection = WENOVectorInvariant(; order=5, time_discretization),
                         tracer_advection = WENO(; order=7, time_discretization, minimum_buffer_upwind_order=1),
                         closure,
                         # forcing = (T = FT, S = FS, u = FTides.u, v = FTides.v),
                         forcing = FTides, 
                         # additional_surface_fluxes = (; T = SurfaceFluxRestoring(FT), S = SurfaceFluxRestoring(FS)), 
                         boundary_conditions = (u = u_obcs, v = v_obcs,
                                                T = T_obcs, S = S_obcs,
                                                # U = U_obcs, V = V_obcs, η = η_obcs, 
                                               ), 
)

# Initialising with GLORYS T and S fields
set!(ocean.model, T = Tᵉˣᵗ[Oceananigans.Units.Time(0.0)], S = Sᵉˣᵗ[Oceananigans.Units.Time(0.0)]) 

# =====================================================================
#  SEA-ICE MODEL (used when POLAR=True)
# =====================================================================
if POLAR
    sea_ice = sea_ice_simulation(grid, ocean; advection = WENO(order=5, minimum_buffer_upwind_order=2))

    ice_dataset = ECCO4Monthly()
    ecco_dir = "../../forcing_data/ecco"
    ecco_h = Metadatum(:sea_ice_thickness; date=dates[1], dataset=ice_dataset, dir=ecco_dir)
    ecco_ℵ = Metadatum(:sea_ice_concentration; date=dates[1], dataset=ice_dataset, dir=ecco_dir)
    set!(sea_ice.model, h=ecco_h, ℵ=ecco_ℵ)
    @info "sea ice ENABLED (initialised from ECCO4)" 
else
    sea_ice = nothing
    @info "sea ice DISABLED (ocean-only)"
end

# =====================================================================
#  ATMOSPHERE (JRA55)
# =====================================================================
atmosphere = JRA55PrescribedAtmosphere(arch; dir_kwJRA55...)
ocean_surface = SurfaceRadiationProperties(albedo = LatitudeDependentAlbedo())
radiation = JRA55PrescribedRadiation(arch; ocean_surface, dir_kwJRA55...)
land = JRA55PrescribedLand(grid; dir_kwJRA55...)

if POLAR
    coupled_model = EarthSystemModel(; ocean, sea_ice, atmosphere, radiation, land)
else
    coupled_model = EarthSystemModel(; ocean, atmosphere, radiation, land)
end

simulation = Simulation(coupled_model; Δt = 5minutes, stop_time = 60days) # initially small time step

wall_time = Ref(time_ns())

function progress(sim)
    ocean = sim.model.ocean
    T = ocean.model.tracers.T
    S = ocean.model.tracers.S
    u, v, w = ocean.model.velocities

    msg = @sprintf("time: %s, iter: %d, extrema(T, S): (%.1f, %.1f) °C (%.1f, %.1f) psu, extrema(u) (%.2e, %.2e, %.2e)",
                   prettytime(sim), iteration(sim),
                   extrema(T)..., extrema(S)..., maximum(abs, u), maximum(abs, v), maximum(abs, w))

    # append sea-ice diagnostics only if ice is in the model
    if POLAR
        h = sim.model.sea_ice.model.ice_thickness
        ℵ = sim.model.sea_ice.model.ice_concentration
        msg *= @sprintf(", max(h): %.2f m, max(ℵ): %.2f", maximum(h), maximum(ℵ))
    end

    msg *= @sprintf(", wall: %s", prettytime(1e-9 * (time_ns() - wall_time[])))

    @info msg
    wall_time[] = time_ns()
    return nothing
end

add_callback!(simulation, progress, TimeInterval(1hours)) # Sets the frequency of the callback

# Saving segment
averaging_schedule = AveragedTimeInterval(6hours, window=6hours)

ocean_outputs = merge(
    ocean.model.tracers,
    ocean.model.velocities,
    (disp = ocean.model.free_surface.displacement, ), 
)

ocean.output_writers[:fields_3d] = JLD2Writer(
    ocean.model, ocean_outputs;
    schedule = averaging_schedule,
    filename = "$(REGION_NAME)_3D_fields_spinup_2010",
    indices = (:, :, :),
    overwrite_files = true,
)

# sea-ice output (only when ice is in the model)
if POLAR
    sea_ice_outputs = (
        h = sea_ice.model.ice_thickness,
        ℵ = sea_ice.model.ice_concentration,
    )
    sea_ice.output_writers[:fields] = JLD2Writer(
        sea_ice.model, sea_ice_outputs;
        schedule = averaging_schedule,
        filename = "$(REGION_NAME)_sea_ice_fields_spinup_2010",
        overwrite_files = true,
    )
    @info "sea-ice output writer added"
end

# turbulent diffusivities
#= function ocean_turbulent_diffusivities(model)
    if hasproperty(model, :closure_fields)
        return model.closure_fields
    elseif hasproperty(model, :diffusivity_fields)
        return model.diffusivity_fields
    else
        return nothing
    end
end

function turbulent_diffusivity_outputs(model)
    cf = ocean_turbulent_diffusivities(model)
    isnothing(cf) && return NamedTuple()

    if cf isa Tuple
        for cfi in cf
            isnothing(cfi) && continue
            if all(Base.Fix1(hasproperty, cfi), (:κu, :κc, :κe))
                return (; κu = cfi.κu, κc = cfi.κc, κe = cfi.κe)
            elseif all(Base.Fix1(hasproperty, cfi), (:κᶜ, :κᵘ))
                return (; κᶜ = cfi.κᶜ, κᵘ = cfi.κᵘ)
            end
        end
        return NamedTuple()
    end

    if all(Base.Fix1(hasproperty, cf), (:κu, :κc, :κe))
        return (; κu = cf.κu, κc = cf.κc, κe = cf.κe)
    elseif all(Base.Fix1(hasproperty, cf), (:κᶜ, :κᵘ))
        return (; κᶜ = cf.κᶜ, κᵘ = cf.κᵘ)
    end

    return NamedTuple()
end

diff_outputs = turbulent_diffusivity_outputs(ocean.model)

if isempty(diff_outputs)
    @warn "No turbulent diffusivity fields found (check closure / Oceananigans version)."
else
    ocean.output_writers[:diffusivities] = JLD2Writer(
        ocean.model, diff_outputs;
        schedule      = averaging_schedule,
        filename      = "$(REGION_NAME)_diffusivity_fields_spinup_2010",
        indices       = (:, :, :),
        overwrite_files = true,
    )
end=#

# Checkpointing for pickups if needed
checkpoint_dir = "$(REGION_NAME)_checkpoints_spinup"
mkpath(checkpoint_dir)

checkpoint_schedule = TimeInterval(7days)

simulation.output_writers[:checkpointer] = Checkpointer(
    simulation.model;
    schedule           = checkpoint_schedule,
    dir                = checkpoint_dir,
    prefix             = "$(REGION_NAME)_sim_ckpt_spinup_2010",
    overwrite_files = true,
    cleanup            = true,
)

# Run simulation
run!(simulation)

using Oceananigans.OutputWriters: write_output!

write_output!(simulation.output_writers[:checkpointer], simulation.model)
