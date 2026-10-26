using Pkg
Pkg.activate(".")
Pkg.instantiate()

DEPOT_PATH

using NumericalEarth, Oceananigans, Oceananigans.Units
using Oceananigans.Units: Time
using Oceananigans.BoundaryConditions: Radiation, FlatherBoundaryCondition, NormalFlowBoundaryCondition
using Oceananigans.Operators: Δzᶠᶜᶜ, Δzᶜᶠᶜ
using Oceananigans.ImmersedBoundaries: immersed_peripheral_node, immersed_inactive_node
using Oceananigans.OrthogonalSphericalShellGrids: RotatedLatitudeLongitudeGrid
using Dates, Printf, Base64
using CopernicusMarine  
using GlobalOceanBioME
using OceanBioME
using Oceananigans.BoundaryConditions: BoundaryConditionOperation

include("hacks.jl")

# Choose the architecture
using CUDA
arch = GPU()
# arch = CPU()

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
const RES   = 1/6        # degrees per cell 1 = 1degree resolution, 1/4 = quarter degree resolution etc.
const DEPTH = 5000meters      # ≥ deepest point in your region
  
# Forcing time window
const DATE_START = DateTime(2010, 9, 1)
const DATE_STOP  = DateTime(2010, 11, 1)
 
# a human label for plot titles / filenames
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

Nz = length(z) # Number of vertical levels
 
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
        halo = (7, 7, 7))
end
 
@info "grid $(Nx)×$(Ny)×$(Nz)  ($(POLAR ? "rotated" : "plain") lat-lon)"

########### # fill in the below with your credentials for NASA Earth Data and Copernicus Marine Data Store

ENV["ECCO_USERNAME"] = 
ENV["ECCO_WEBDAV_PASSWORD"] = 
ENV["COPERNICUS_USERNAME"] = 
ENV["COPERNICUS_PASSWORD"] = 

########### Set these to the directories you want the forcing data to download to: GLORYS is the boundary 
########### forcing and JRA55 is the atmospheric forcing
########### This also now includes a directory for biological files from ECCO Darwin

ENV["GLORYS_DIR"] = "Forcings/GLORYSDaily"
ENV["JRA55_DIR"] = "Forcings/JRA55"
ENV["ECCODarwin_DIR"] = "Forcings/ECCODarwin"

dir_kw = haskey(ENV, "DATA_DIR")   ? (; dir = ENV["DATA_DIR"])   : (;)
dir_kwGLORYS = haskey(ENV, "GLORYS_DIR") ? (; dir = ENV["GLORYS_DIR"]) : (;)
dir_kwJRA55 = haskey(ENV, "JRA55_DIR")  ? (; dir = ENV["JRA55_DIR"])  : (;)
dir_kwECCODarwin = haskey(ENV, "ECCODarwin_DIR")  ? (; dir = ENV["ECCODarwin_DIR"])  : (;)

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
using NCDatasets, Interpolations

ds = Dataset("Topography/topography.nc") # Replace with your downloaded file name
raw_lon = ds["lon"][:]
raw_lat = ds["lat"][:]
raw_elevation = ds["elevation"][:, :]

λ = λnodes(underlying_grid, Center())
φ = φnodes(underlying_grid, Center())

gebco_interp = LinearInterpolation((raw_lon, raw_lat), raw_elevation)
bathymetry_matrix = [gebco_interp(l, p) for l in λ, p in φ]

grid = ImmersedBoundaryGrid(underlying_grid, PartialCellBottom(bathymetry_matrix); active_cells_map = true)

#  Plot the bathymetry / land mask to check the region.
using CairoMakie, SixelTerm
 
h_bottom = Array(interior(grid.immersed_boundary.bottom_height, :, :, 1))
h_bottom[h_bottom .≥ 0] .= NaN
 
fig = Figure(size = (700, 600))
ax = Axis(fig[1, 1],
          xlabel = POLAR ? "grid i (rotated)" : "longitude index",
          ylabel = POLAR ? "grid j (rotated)" : "latitude index",
          title  = "$REGION_NAME bathymetry ($(POLAR ? "rotated" : "plain") grid)",
          aspect = DataAspect())
hm = heatmap!(ax, h_bottom, colormap = :deep, colorrange = (-DEPTH, 0), nan_color = :gray80)
Colorbar(fig[1, 2], hm, label = "bottom height [m]")
save("bathymetry.pdf", fig, pdf_version="1.4")


# =====================================================================
#  ATMOSPHERE 
# =====================================================================
# JRA55_kw = (; dataset = MultiYearJRA55(), region = BoundingBox(longitude = (0, 360), latitude = (φ₁ - 5, φ₂ + 5)), start_date = DATE_START, end_date = DATE_STOP, dir_kwJRA55...)
JRA55_kw = (; region = BoundingBox(longitude = (0, 360), latitude = (φ₁ - 5, φ₂ + 5)), dir_kwJRA55...)

atmosphere = JRA55PrescribedAtmosphere(arch; JRA55_kw...)
ocean_surface = SurfaceRadiationProperties(albedo = LatitudeDependentAlbedo())
radiation = JRA55PrescribedRadiation(arch; ocean_surface, JRA55_kw...)
land = JRA55PrescribedLand(arch; JRA55_kw...)

# Which edges are OPEN ocean boundaries?
# The open boundaries will be forced with the GLORYSDaily dataset (daily temperature, salinity and velocities)
# Or from the ECCO4-Darwin dataset (for the Biogeochemical fields)
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

# ---- forcing time series (GLORYS), region tracks the domain ---------
dates   = DATE_START : Day(1) : DATE_STOP
dataset = GLORYSDaily()
# pad the forcing region slightly beyond the model bounds so edges are covered
region  = BoundingBox(longitude = (λ₁ - 1, λ₂ + 1), latitude = (φ₁ - 1, φ₂ + 1))

Tᵉˣᵗ = FieldTimeSeries(Metadata(:temperature; dates, dataset, region, dir_kwGLORYS...), grid, inpainting=100, time_indices_in_memory=2)
Sᵉˣᵗ = FieldTimeSeries(Metadata(:salinity; dates, dataset, region, dir_kwGLORYS...), grid, inpainting=100, time_indices_in_memory=2)
uᵉˣᵗ = FieldTimeSeries(Metadata(:u_velocity; dates, dataset, region, dir_kwGLORYS...), grid, inpainting=100, time_indices_in_memory=2)
vᵉˣᵗ = FieldTimeSeries(Metadata(:v_velocity; dates, dataset, region, dir_kwGLORYS...), grid, inpainting=100, time_indices_in_memory=2)
ηᵉˣᵗ = FieldTimeSeries(Metadata(:free_surface; dates, dataset, region, dir_kwGLORYS...), grid, inpainting=100, time_indices_in_memory=2)
nothing #hide

# ---- BGC forcing dataset (ECCO4-Darwin, monthly) --------------------
bgc_dataset = ECCO4DarwinMonthly()
bgc_dates   = DATE_START : Month(1) : DATE_STOP
bgc_kwargs  = (dataset = bgc_dataset, dates = bgc_dates, region = region)

# ---- external BGC fields (all present in ECCO4-Darwin) --------------
NO₃ᵉˣᵗ = FieldTimeSeries(Metadata(:nitrate; bgc_kwargs..., dir_kwECCODarwin...), grid, inpainting=100, time_indices_in_memory=2)
NH₄ᵉˣᵗ = FieldTimeSeries(Metadata(:ammonia; bgc_kwargs..., dir_kwECCODarwin...), grid, inpainting=100, time_indices_in_memory=2)
PO₄ᵉˣᵗ = FieldTimeSeries(Metadata(:phosphate; bgc_kwargs..., dir_kwECCODarwin...), grid, inpainting=100, time_indices_in_memory=2)
Siᵉˣᵗ  = FieldTimeSeries(Metadata(:silicate; bgc_kwargs..., dir_kwECCODarwin...), grid, inpainting=100, time_indices_in_memory=2)
DICᵉˣᵗ = FieldTimeSeries(Metadata(:dissolved_inorganic_carbon; bgc_kwargs..., dir_kwECCODarwin...), grid, inpainting=100, time_indices_in_memory=2)
Alkᵉˣᵗ = FieldTimeSeries(Metadata(:alk; bgc_kwargs..., dir_kwECCODarwin...), grid, inpainting=100, time_indices_in_memory=2)
O2ᵉˣᵗ = FieldTimeSeries(Metadata(:dissolved_oxygen; bgc_kwargs..., dir_kwECCODarwin...), grid, inpainting=100, time_indices_in_memory=2)
POCᵉˣᵗ = FieldTimeSeries(Metadata(:particulate_organic_carbon; bgc_kwargs..., dir_kwECCODarwin...), grid, inpainting=100, time_indices_in_memory=2)
DOCᵉˣᵗ = FieldTimeSeries(Metadata(:dissolved_organic_carbon; bgc_kwargs..., dir_kwECCODarwin...), grid, inpainting=100, time_indices_in_memory=2)
Pᵉˣᵗ = FieldTimeSeries(Metadata(:synechococcus; bgc_kwargs..., dir_kwECCODarwin...), grid, inpainting=100, time_indices_in_memory=2)
Zᵉˣᵗ = FieldTimeSeries(Metadata(:zooplankton_grazing_small; bgc_kwargs..., dir_kwECCODarwin...), grid, inpainting=100, time_indices_in_memory=2)
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

# ---- tracer boundary functions (one per edge) -----------------------
# west/east are i-edges → the free index is j; south/north are j-edges → free index is i.
@inline west_trc( j, k, grid, clock, fields, φ) = @inbounds φ[1,       j, k, Oceananigans.Units.Time(clock.time)]
@inline east_trc( j, k, grid, clock, fields, φ) = @inbounds φ[grid.Nx, j, k, Oceananigans.Units.Time(clock.time)]
@inline south_trc(i, k, grid, clock, fields, φ) = @inbounds φ[i, 1,       k, Oceananigans.Units.Time(clock.time)]
@inline north_trc(i, k, grid, clock, fields, φ) = @inbounds φ[i, grid.Ny, k, Oceananigans.Units.Time(clock.time)]

# ---- normal-velocity boundary functions -----------------------------
# west/east: normal component is u, indexed at the i-face
@inline west_u( j, k, grid, clock, fields, φ) = @inbounds φ[1,        j, k, Oceananigans.Units.Time(clock.time)]
@inline east_u( j, k, grid, clock, fields, φ) = @inbounds φ[grid.Nx+1, j, k, Oceananigans.Units.Time(clock.time)]
# south/north: normal component is v, indexed at the j-face
@inline south_v(i, k, grid, clock, fields, φ) = @inbounds φ[i, 1,        k, Oceananigans.Units.Time(clock.time)]
@inline north_v(i, k, grid, clock, fields, φ) = @inbounds φ[i, grid.Ny+1, k, Oceananigans.Units.Time(clock.time)]

# ---- Flather barotropic functions (return (transport, η)) -----------
@inline function west_U(j, k, grid, clock, fields, p)
    t = Oceananigans.Units.Time(clock.time)
    U = vertical_integral(1, j, grid, p.u, t, Δzᶠᶜᶜ, Face(), Center(), Center())
    return (U, @inbounds p.η[1, j, 1, t])
end
@inline function east_U(j, k, grid, clock, fields, p)
    t = Oceananigans.Units.Time(clock.time); i = grid.Nx
    U = vertical_integral(i, j, grid, p.u, t, Δzᶠᶜᶜ, Face(), Center(), Center())
    return (U, @inbounds p.η[i, j, 1, t])
end
@inline function south_V(i, k, grid, clock, fields, p)
    t = Oceananigans.Units.Time(clock.time)
    V = vertical_integral(i, 1, grid, p.v, t, Δzᶜᶠᶜ, Center(), Face(), Center())
    return (V, @inbounds p.η[i, 1, 1, t])
end
@inline function north_V(i, k, grid, clock, fields, p)
    t = Oceananigans.Units.Time(clock.time); j = grid.Ny
    V = vertical_integral(i, j, grid, p.v, t, Δzᶜᶠᶜ, Center(), Face(), Center())
    return (V, @inbounds p.η[i, j, 1, t])
end

# ---- assemble OBC keyword pairs only for the OPEN edges -------------
rad = Radiation(inflow_timescale = 1days)

T_pairs = Pair{Symbol,Any}[]; S_pairs = Pair{Symbol,Any}[]
u_pairs = Pair{Symbol,Any}[]; v_pairs = Pair{Symbol,Any}[]
U_pairs = Pair{Symbol,Any}[]; V_pairs = Pair{Symbol,Any}[]
# BGC tracer OBC pair-lists (data-fed from Darwin on the open edges)
NO3_pairs = Pair{Symbol,Any}[]; NH4_pairs = Pair{Symbol,Any}[]; P_pairs = Pair{Symbol,Any}[]; Si_pairs = Pair{Symbol,Any}[]
DIC_pairs = Pair{Symbol,Any}[]; Alk_pairs = Pair{Symbol,Any}[]; O2_pairs = Pair{Symbol,Any}[]
POC_pairs = Pair{Symbol,Any}[]; DOC_pairs = Pair{Symbol,Any}[]
P_pairs = Pair{Symbol,Any}[]; Z_pairs = Pair{Symbol,Any}[]

if get(OPEN, :west, false)
    push!(T_pairs, :west => ValueBoundaryCondition(west_trc, discrete_form=true, parameters=Tᵉˣᵗ, scheme=rad))
    push!(S_pairs, :west => ValueBoundaryCondition(west_trc, discrete_form=true, parameters=Sᵉˣᵗ, scheme=rad))
    push!(NO3_pairs,   :west => ValueBoundaryCondition(west_trc, discrete_form=true, parameters=NO₃ᵉˣᵗ, scheme=rad))
    push!(NH4_pairs,   :west => ValueBoundaryCondition(west_trc, discrete_form=true, parameters=NH₄ᵉˣᵗ, scheme=rad))
    push!(P_pairs,   :west => ValueBoundaryCondition(west_trc, discrete_form=true, parameters=PO₄ᵉˣᵗ, scheme=rad))
    push!(Si_pairs,  :west => ValueBoundaryCondition(west_trc, discrete_form=true, parameters=Siᵉˣᵗ, scheme=rad))
    push!(DIC_pairs, :west => ValueBoundaryCondition(west_trc, discrete_form=true, parameters=DICᵉˣᵗ, scheme=rad))
    push!(Alk_pairs, :west => ValueBoundaryCondition(west_trc, discrete_form=true, parameters=Alkᵉˣᵗ, scheme=rad))
    push!(O2_pairs, :west => ValueBoundaryCondition(west_trc, discrete_form=true, parameters=O2ᵉˣᵗ, scheme=rad))
    push!(POC_pairs, :west => ValueBoundaryCondition(west_trc, discrete_form=true, parameters=POCᵉˣᵗ, scheme=rad))
    push!(DOC_pairs, :west => ValueBoundaryCondition(west_trc, discrete_form=true, parameters=DOCᵉˣᵗ, scheme=rad))
    push!(P_pairs, :west => ValueBoundaryCondition(west_trc, discrete_form=true, parameters=Pᵉˣᵗ, scheme=rad))
    push!(Z_pairs, :west => ValueBoundaryCondition(west_trc, discrete_form=true, parameters=Zᵉˣᵗ, scheme=rad))
    push!(u_pairs, :west => NormalFlowBoundaryCondition(west_u, discrete_form=true, parameters=uᵉˣᵗ, scheme=rad))
    push!(U_pairs, :west => FlatherBoundaryCondition(west_U, discrete_form=true, parameters=(u=uᵉˣᵗ, η=ηᵉˣᵗ)))
end
if get(OPEN, :east, false)
    push!(T_pairs, :east => ValueBoundaryCondition(east_trc, discrete_form=true, parameters=Tᵉˣᵗ, scheme=rad))
    push!(S_pairs, :east => ValueBoundaryCondition(east_trc, discrete_form=true, parameters=Sᵉˣᵗ, scheme=rad))
    push!(NO3_pairs,   :east => ValueBoundaryCondition(east_trc, discrete_form=true, parameters=NO₃ᵉˣᵗ, scheme=rad))
    push!(NH4_pairs,   :east => ValueBoundaryCondition(east_trc, discrete_form=true, parameters=NH₄ᵉˣᵗ, scheme=rad))
    push!(P_pairs,   :east => ValueBoundaryCondition(east_trc, discrete_form=true, parameters=PO₄ᵉˣᵗ, scheme=rad))
    push!(Si_pairs,  :east => ValueBoundaryCondition(east_trc, discrete_form=true, parameters=Siᵉˣᵗ, scheme=rad))
    push!(DIC_pairs, :east => ValueBoundaryCondition(east_trc, discrete_form=true, parameters=DICᵉˣᵗ, scheme=rad))
    push!(Alk_pairs, :east => ValueBoundaryCondition(east_trc, discrete_form=true, parameters=Alkᵉˣᵗ, scheme=rad))
    push!(O2_pairs, :east => ValueBoundaryCondition(east_trc, discrete_form=true, parameters=O2ᵉˣᵗ, scheme=rad))
    push!(POC_pairs, :east => ValueBoundaryCondition(east_trc, discrete_form=true, parameters=POCᵉˣᵗ, scheme=rad))
    push!(DOC_pairs, :east => ValueBoundaryCondition(east_trc, discrete_form=true, parameters=DOCᵉˣᵗ, scheme=rad))
    push!(P_pairs, :east => ValueBoundaryCondition(east_trc, discrete_form=true, parameters=Pᵉˣᵗ, scheme=rad))
    push!(Z_pairs, :east => ValueBoundaryCondition(east_trc, discrete_form=true, parameters=Zᵉˣᵗ, scheme=rad))
    push!(u_pairs, :east => NormalFlowBoundaryCondition(east_u, discrete_form=true, parameters=uᵉˣᵗ, scheme=rad))
    push!(U_pairs, :east => FlatherBoundaryCondition(east_U, discrete_form=true, parameters=(u=uᵉˣᵗ, η=ηᵉˣᵗ)))
end
if get(OPEN, :south, false)
    push!(T_pairs, :south => ValueBoundaryCondition(south_trc, discrete_form=true, parameters=Tᵉˣᵗ, scheme=rad))
    push!(S_pairs, :south => ValueBoundaryCondition(south_trc, discrete_form=true, parameters=Sᵉˣᵗ, scheme=rad))
    push!(NO3_pairs,   :south => ValueBoundaryCondition(south_trc, discrete_form=true, parameters=NO₃ᵉˣᵗ, scheme=rad))
    push!(NH4_pairs,   :south => ValueBoundaryCondition(south_trc, discrete_form=true, parameters=NH₄ᵉˣᵗ, scheme=rad))
    push!(P_pairs,   :south => ValueBoundaryCondition(south_trc, discrete_form=true, parameters=PO₄ᵉˣᵗ, scheme=rad))
    push!(Si_pairs,  :south => ValueBoundaryCondition(south_trc, discrete_form=true, parameters=Siᵉˣᵗ, scheme=rad))
    push!(DIC_pairs, :south => ValueBoundaryCondition(south_trc, discrete_form=true, parameters=DICᵉˣᵗ, scheme=rad))
    push!(Alk_pairs, :south => ValueBoundaryCondition(south_trc, discrete_form=true, parameters=Alkᵉˣᵗ, scheme=rad))
    push!(O2_pairs, :south => ValueBoundaryCondition(south_trc, discrete_form=true, parameters=O2ᵉˣᵗ, scheme=rad))
    push!(POC_pairs, :south => ValueBoundaryCondition(south_trc, discrete_form=true, parameters=POCᵉˣᵗ, scheme=rad))
    push!(DOC_pairs, :south => ValueBoundaryCondition(south_trc, discrete_form=true, parameters=DOCᵉˣᵗ, scheme=rad))
    push!(P_pairs, :south => ValueBoundaryCondition(south_trc, discrete_form=true, parameters=Pᵉˣᵗ, scheme=rad))
    push!(Z_pairs, :south => ValueBoundaryCondition(south_trc, discrete_form=true, parameters=Zᵉˣᵗ, scheme=rad))
    push!(v_pairs, :south => NormalFlowBoundaryCondition(south_v, discrete_form=true, parameters=vᵉˣᵗ, scheme=rad))
    push!(V_pairs, :south => FlatherBoundaryCondition(south_V, discrete_form=true, parameters=(v=vᵉˣᵗ, η=ηᵉˣᵗ)))
end
if get(OPEN, :north, false)
    push!(T_pairs, :north => ValueBoundaryCondition(north_trc, discrete_form=true, parameters=Tᵉˣᵗ, scheme=rad))
    push!(S_pairs, :north => ValueBoundaryCondition(north_trc, discrete_form=true, parameters=Sᵉˣᵗ, scheme=rad))
    push!(NO3_pairs,   :north => ValueBoundaryCondition(north_trc, discrete_form=true, parameters=NO₃ᵉˣᵗ, scheme=rad))
    push!(NH4_pairs,   :north => ValueBoundaryCondition(north_trc, discrete_form=true, parameters=NH₄ᵉˣᵗ, scheme=rad))
    push!(P_pairs,   :north => ValueBoundaryCondition(north_trc, discrete_form=true, parameters=PO₄ᵉˣᵗ, scheme=rad))
    push!(Si_pairs,  :north => ValueBoundaryCondition(north_trc, discrete_form=true, parameters=Siᵉˣᵗ, scheme=rad))
    push!(DIC_pairs, :north => ValueBoundaryCondition(north_trc, discrete_form=true, parameters=DICᵉˣᵗ, scheme=rad))
    push!(Alk_pairs, :north => ValueBoundaryCondition(north_trc, discrete_form=true, parameters=Alkᵉˣᵗ, scheme=rad))
    push!(O2_pairs, :north => ValueBoundaryCondition(north_trc, discrete_form=true, parameters=O2ᵉˣᵗ, scheme=rad))
    push!(POC_pairs, :north => ValueBoundaryCondition(north_trc, discrete_form=true, parameters=POCᵉˣᵗ, scheme=rad))
    push!(DOC_pairs, :north => ValueBoundaryCondition(north_trc, discrete_form=true, parameters=DOCᵉˣᵗ, scheme=rad))
    push!(P_pairs, :north => ValueBoundaryCondition(north_trc, discrete_form=true, parameters=Pᵉˣᵗ, scheme=rad))
    push!(Z_pairs, :north => ValueBoundaryCondition(north_trc, discrete_form=true, parameters=Zᵉˣᵗ, scheme=rad))
    push!(v_pairs, :north => NormalFlowBoundaryCondition(north_v, discrete_form=true, parameters=vᵉˣᵗ, scheme=rad))
    push!(V_pairs, :north => FlatherBoundaryCondition(north_V, discrete_form=true, parameters=(v=vᵉˣᵗ, η=ηᵉˣᵗ)))
end

T_obcs = FieldBoundaryConditions(; T_pairs...)
S_obcs = FieldBoundaryConditions(; S_pairs...)
u_obcs = FieldBoundaryConditions(; u_pairs...)
v_obcs = FieldBoundaryConditions(; v_pairs...)
U_obcs = FieldBoundaryConditions(grid, (Face(),   Center(), nothing); U_pairs...)
V_obcs = FieldBoundaryConditions(grid, (Center(), Face(),   nothing); V_pairs...)

# ---- BGC tracer OBCs ------------------------------------------------
NO3_obcs   = FieldBoundaryConditions(; NO3_pairs...)
NH4_obcs   = FieldBoundaryConditions(; NH4_pairs...)
PO₄_obcs = FieldBoundaryConditions(; P_pairs...)
Si_obcs  = FieldBoundaryConditions(; Si_pairs...)
Alk_obcs = FieldBoundaryConditions(; Alk_pairs...)
POC_obcs = FieldBoundaryConditions(; POC_pairs...)
DOC_obcs = FieldBoundaryConditions(; DOC_pairs...)
P_obcs = FieldBoundaryConditions(; P_pairs...)
Z_obcs = FieldBoundaryConditions(; Z_pairs...)

# =====================================================================
#  SPONGE LAYERS 
# =====================================================================
# We also add a sponge layers at open edges to damp any spurious waves. 
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

sponge_width_TS = 2   # cells
sponge_mask_TS  = zeros(Nx, Ny, Nz)
for j in 1:Ny, i in 1:Nx, k in 1:Nz
    dists = Float64[]
    get(OPEN, :west,  false) && push!(dists, i - 1)
    get(OPEN, :east,  false) && push!(dists, Nx - i)
    get(OPEN, :south, false) && push!(dists, j - 1)
    get(OPEN, :north, false) && push!(dists, Ny - j)
    push!(dists, Nz - k)
    d = isempty(dists) ? Inf : minimum(dists)      # nearest open edge
    sponge_mask_TS[i, j, k] = exp(-(d / sponge_width_TS)^2)
end

FT = DatasetRestoring(Metadata(:temperature; dates, dataset, region, dir_kwGLORYS...), grid; rate = 1/1days, mask = sponge_mask_TS, inpainting=100)
Fu = DatasetRestoring(Metadata(:u_velocity; dates, dataset, region, dir_kwGLORYS...), grid; rate = 1/20minutes, mask = sponge_mask, inpainting=100)
Fv = DatasetRestoring(Metadata(:v_velocity; dates, dataset, region, dir_kwGLORYS...), grid; rate = 1/20minutes, mask = sponge_mask, inpainting=100)
FS = DatasetRestoring(Metadata(:salinity; dates, dataset, region, dir_kwGLORYS...), grid; rate = 1/1days, mask = sponge_mask_TS, inpainting=100)

# ---- BGC sponge restoring (same mask) ------------------------------
FNO₃ = DatasetRestoring(Metadata(:nitrate; bgc_kwargs..., dir_kwECCODarwin...), grid; rate = 1/1days, mask = sponge_mask, inpainting=100)
FNH₄ = DatasetRestoring(Metadata(:ammonia; bgc_kwargs..., dir_kwECCODarwin...), grid; rate = 1/1days, mask = sponge_mask, inpainting=100)
FSi  = DatasetRestoring(Metadata(:silicate; bgc_kwargs..., dir_kwECCODarwin...), grid; rate = 1/1days, mask = sponge_mask, inpainting=100)
FPO₄ = DatasetRestoring(Metadata(:phosphate; bgc_kwargs..., dir_kwECCODarwin...), grid; rate = 1/1days, mask = sponge_mask, inpainting=100)
FDIC = DatasetRestoring(Metadata(:dissolved_inorganic_carbon; bgc_kwargs..., dir_kwECCODarwin...), grid; rate = 1/1days, mask = sponge_mask, inpainting=100)
FAlk = DatasetRestoring(Metadata(:alk; bgc_kwargs..., dir_kwECCODarwin...), grid; rate = 1/1days, mask = sponge_mask, inpainting=100)
FO2 = DatasetRestoring(Metadata(:dissolved_oxygen; bgc_kwargs..., dir_kwECCODarwin...), grid; rate = 1/1days, mask = sponge_mask, inpainting=100)
FPOC = DatasetRestoring(Metadata(:particulate_organic_carbon; bgc_kwargs..., dir_kwECCODarwin...), grid; rate = 1/1days, mask = sponge_mask, inpainting=100)
FDOC = DatasetRestoring(Metadata(:dissolved_organic_carbon; bgc_kwargs..., dir_kwECCODarwin...), grid; rate = 1/1days, mask = sponge_mask, inpainting=100)
FP = DatasetRestoring(Metadata(:synechococcus; bgc_kwargs..., dir_kwECCODarwin...), grid; rate = 1/1days, mask = sponge_mask, inpainting=100)
FZ = DatasetRestoring(Metadata(:zooplankton_grazing_small; bgc_kwargs..., dir_kwECCODarwin...), grid; rate = 1/1days, mask = sponge_mask, inpainting=100)

@info "open edges: $(join([string(e) for (e,v) in pairs(OPEN) if v], ", "))"

# Mixing closure
closure = CATKEVerticalDiffusivity(minimum_tke=1e-7)
time_discretization = AdaptiveVerticallyImplicitDiscretization(cfl=0.5)

# =====================================================================
#  BIOGEOCHEMISTRY
# =====================================================================

# Light for biogeochemistry
if POLAR
    underlying_light_attenuation = fixed_attenuation_par_from_radiation(grid, radiation; PAR_fraction = 0.43, attenuation_coefficient = 0.3)
    light_attenuation = GlobalOceanBioME.IceMaskedLightAttenuation(; underlying_light_attenuation,
                                                                     ice_thickness = Field{Center, Center, Nothing}(grid),
                                                                     ice_concentration = Field{Center, Center, Nothing}(grid))
else
    # using Oceananigans.Units: Time
    
    # light_attenuation = fixed_attenuation_par_from_radiation(grid, radiation; PAR_fraction = 1.0, attenuation_coefficient = 0.01)

    @inline function surface_PAR(x, y, t)
        pos_sun_360 = 360 / (24hours) * t % 360
        pos_sun = (pos_sun_360 + 180) % 360 - 180
        return 50 * exp(-(x - pos_sun)^2 / 500)
    end
    
    light_attenuation = PrescribedAttenuationPAR(grid, surface_PAR; attenuation = 0.04)
end

# Conservation laws for Plankton breaking down organic matter into inorganic one, 
# Nitrients consumed by P and Z that then become organic matter, 
# Remineralisation of organic matter into inorganic one in the presence of Oxygen (oxidant) 
# and Nitrification in the presence of Oxygen, respectively
modifiers = (ScaleNegativeTracers((:P, :Z, :DOM, :POM, :DIC); invalid_fill_value = zero(grid)), 
             ScaleNegativeTracers((:NO₃, :NH₄, :PO₄, :Si, :P, :Z, :DOM, :POM); invalid_fill_value = zero(grid)), 
             ScaleNegativeTracers((:O₂, :DOM, :POM, :DIC); invalid_fill_value = zero(grid)), 
             ScaleNegativeTracers((:O₂, :NO₃, :NH₄); invalid_fill_value = zero(grid))) 

# Detritus model
detritus = DissolvedParticulate(grid, :DOM, :POM;
                              dissolved_remineralisation_rate = 1 / (6 * 30 * days), # 1/s
                              particulate_remineralisation_rate = 0.03 / days,
                              dissolved_fraction_of_remineralisation = 0.0,
                              sinking_speeds = 10 / days,
                              open_bottom = true)

# Plankton model
using Oceananigans.Fields: ZeroField
# limiting_nutrients = (:nitrate, :ammonia)
plankton = PhytoZoo(grid; edible_detritus_name = :POM, # zooplankton_assimilation_fraction = 1.0, 
                  temperature_coefficient = 1.88, 
                  phytoplankton_mortality_rate = 1e-6,
                  # zooplankton_mortality_rate = 1e-7, 
                  # nutrient_half_saturations = (nitrate = 2.3, ammonia = 0.001)[limiting_nutrients], 
                  phytoplankton_sinking_velocity = (u = ZeroField(), v = ZeroField(), w = ZeroField()),
                  zooplankton_sinking_velocity = (u = ZeroField(), v = ZeroField(), w = ZeroField()))

# Biogeochemistry model that uses the NutrientPlanktonDetritus framework of OceanBioME.jl
biogeochemistry = NutrientsPlanktonDetritus(grid; light_attenuation, nutrients = Nutrients(NitrateAmmonia(), PO₄, nothing, Si), plankton = plankton, detritus = detritus, inorganic_carbon = CarbonateSystem(), oxygen = Oxygen(), modifiers)

@info Oceananigans.Biogeochemistry.required_biogeochemical_tracers(biogeochemistry)

# air–sea CO₂ and O₂ exchange fluxes
underlying_CO₂_flux = CarbonDioxideGasExchangeBoundaryCondition(; wind_speed = 2, 
                                                                air_concentration = 357.21)
# CO₂_flux = GlobalOceanBioME.IceMaskedGasExchangeBoundaryCondition(underlying_CO₂_flux, Field{Center, Center, Nothing}(grid))
push!(DIC_pairs, :top => underlying_CO₂_flux)
underlying_O₂_flux = OxygenGasExchangeBoundaryCondition(; wind_speed = 2)
push!(O2_pairs, :top => underlying_O₂_flux)
DIC_obcs = FieldBoundaryConditions(; DIC_pairs...)
O2_obcs = FieldBoundaryConditions(; O2_pairs...)

# @info DIC_obcs

# =====================================================================
#  OCEAN MODEL (WITH BIOGEOCHEMISTRY)
# =====================================================================
ocean = ocean_simulation(grid;
                         biogeochemistry,
                         free_surface = SplitExplicitFreeSurface(grid; substeps=80),
                         momentum_advection = WENOVectorInvariant(; order=5, time_discretization),
			             tracer_advection = WENO(; order = 5, time_discretization),
                         closure,
                         tracers = (:NO₃, :NH₄, :PO₄, :Si, :P, :Z, :T, :S, :DOM, :POM, :DIC, :Alk, :O₂, ),
                         forcing = (T = FT, S = FS, u = Fu, v = Fv, 
                                    NO₃ = FNO₃, NH₄ = FNH₄, PO₄ = FPO₄, Si = FSi, 
                                    DIC = FDIC, Alk = FAlk, O₂ = FO2, 
                                    POM = FPOC, DOM = FDOC, 
                                    P = FP, Z = FZ, 
                                    ),
                         boundary_conditions = (u = u_obcs, v = v_obcs, T = T_obcs, S = S_obcs, U = U_obcs, V = V_obcs,
                                                NO₃ = NO3_obcs, NH₄ = NH4_obcs, PO₄ = PO₄_obcs, Si = Si_obcs, 
                                                DIC = DIC_obcs, Alk = Alk_obcs, O₂ = O2_obcs,
                                                POM = POC_obcs, DOM = DOC_obcs,
                                                P = P_obcs, Z = Z_obcs, 
                                                )
)     

# initial conditions: physics from GLORYS, BGC from ECCO4-Darwin
set!(ocean.model, T = Tᵉˣᵗ[1], S = Sᵉˣᵗ[1],
                  NO₃ = NO₃ᵉˣᵗ[1], NH₄ = NH₄ᵉˣᵗ[1], PO₄ = PO₄ᵉˣᵗ[1], Si = Siᵉˣᵗ[1], 
                  DIC = DICᵉˣᵗ[1], Alk = Alkᵉˣᵗ[1], O₂ = O2ᵉˣᵗ[1], 
                  POM = POCᵉˣᵗ[1], DOM = DOCᵉˣᵗ[1], 
                  P = Pᵉˣᵗ[1], Z = Zᵉˣᵗ[1])

# =====================================================================
#  SEA-ICE MODEL (used when POLAR=True)
# =====================================================================

# sea ice (only when POLAR)
if POLAR
    sea_ice = sea_ice_simulation(grid, ocean; advection = WENO(order=5, minimum_buffer_upwind_order=2))

    ice_dataset = ECCO4Monthly()
    ecco_dir = "../../forcing_data/ecco"
    ecco_h = Metadatum(:sea_ice_thickness;     date=dates[1], dataset=ice_dataset, dir=ecco_dir)
    ecco_ℵ = Metadatum(:sea_ice_concentration; date=dates[1], dataset=ice_dataset, dir=ecco_dir)
    set!(sea_ice.model, h=ecco_h, ℵ=ecco_ℵ)
    @info "sea ice ENABLED (initialised from ECCO4)"
else
    sea_ice = nothing
    @info "sea ice DISABLED (ocean-only)"
end

# =====================================================================
#  COUPLE MODEL
# =====================================================================
if POLAR
    coupled_model = EarthSystemModel(; ocean, sea_ice, land, atmosphere, radiation)
else
    coupled_model = EarthSystemModel(; ocean, land, atmosphere, radiation)
end

# Change the Δt and simulation length here
simulation = Simulation(coupled_model; Δt = 5minutes, stop_time = 60days) # initially small time step

# This sets the output you see during the simulation
wall_time = Ref(time_ns())

function progress(sim)
    ocean = sim.model.ocean
    T = ocean.model.tracers.T
    S = ocean.model.tracers.S
    u, v, w = ocean.model.velocities

    msg = @sprintf("time: %s, iter: %d, extrema(T, S): (%.1f, %.1f) °C (%.1f, %.1f) psu, extrema(u) (%.2e, %.2e, %.2e)",
                   prettytime(sim), iteration(sim),
                   extrema(T)..., extrema(S)..., maximum(abs, u), maximum(abs, v), maximum(abs, w))

    # biogeochemical tracers
    # Printed only if present, so this cell still works if you change
    # the biogeochemistry or drop a tracer.
    tracers = ocean.model.tracers
    for (name, label, fmt) in (("Z",   "Z",   "%.2f"),
                               ("P",   "P",   "%.2f"),
                               ("NO₃",   "NO₃",   "%.2f"),
                               ("NH₄",   "NH₄",   "%.2f"),
                               ("DIC", "DIC", "%.0f"),
                               ("Alk", "Alk", "%.0f"), 
                               ("O₂", "O₂", "%.0f"), 
                               ("Si", "Si", "%.2f"),
                               ("DOM", "DOM", "%.2f"),
                               ("sPOM", "sPOM", "%.2f"),
                               ("bPOM", "bPOM", "%.2f"),
                               ("POM", "POM", "%.2f"),
                               ("PO₄", "PO₄", "%.2f"))
        sym = Symbol(name)
        if hasproperty(tracers, sym)
            lo, hi = extrema(getproperty(tracers, sym))
            msg *= ", " * label * @eval(@sprintf($(" ∈ [" * fmt * ", " * fmt * "]"), $lo, $hi))
        end
    end

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

add_callback!(simulation, progress, TimeInterval(4hours))

# Saving segment
averaging_schedule = AveragedTimeInterval(6hours, window=6hours)

# air–sea CO₂ flux diagnostic (top BC operation on DIC)
CO₂_flux = Oceananigans.Models.BoundaryConditionOperation(ocean.model.tracers.DIC, :top, ocean.model)

ocean_outputs = merge(
    ocean.model.tracers,          # now includes N, PO₄, Fe, DIC, Alk (+ T,S,e)
    ocean.model.velocities,
)

ocean.output_writers[:fields_3d] = JLD2Writer(
    ocean.model, ocean_outputs;
    schedule = averaging_schedule,
    filename = "$(REGION_NAME)_3D_fields_spinup_$(Dates.year(DATE_START))",
    indices = (:, :, :),
    overwrite_existing = true,
)

#=ocean.output_writers[:PAR] = JLD2Writer(
    ocean.model, (; PAR = ocean.model.biogeochemistry.light_attenuation.field);
    schedule = TimeInterval(1hours), 
    filename = "$(REGION_NAME)_PAR_spinup_$(Dates.year(DATE_START))",
    indices = (:, :, :),
    overwrite_existing = true,
)=#

# turbulent diffusivities 
function ocean_turbulent_diffusivities(model)
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
        filename      = "$(REGION_NAME)_diffusivity_fields_spinup_$(Dates.year(DATE_START))",
        indices       = (:, :, :),
        overwrite_existing = true,
    )
end

# Checkpointing for pickups if needed
checkpoint_dir = "$(REGION_NAME)_checkpoints_spinup"
mkpath(checkpoint_dir)

checkpoint_schedule = TimeInterval(1days)

simulation.output_writers[:checkpointer] = Checkpointer(
    simulation.model;
    schedule           = checkpoint_schedule,
    dir                = checkpoint_dir,
    prefix             = "$(REGION_NAME)_sim_ckpt_spinup_$(Dates.year(DATE_START))",
    overwrite_existing = true,
    cleanup            = true,
)

# Run the simulation
run!(simulation)

using Oceananigans.OutputWriters: write_output!

write_output!(simulation.output_writers[:checkpointer], simulation.model)

