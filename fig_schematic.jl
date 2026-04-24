using CairoMakie
using LaTeXStrings

include(joinpath(@__DIR__, "fig_schematic_geometry.jl"))

# =========================================================
# Schematic: two-layer mode-conversion geometry
# Self-contained — no external dependencies
# =========================================================

fig = Figure(size=(520, 420), fontsize=18)

ax = Axis(fig[1,1],
    aspect = DataAspect(),
    limits = (-2.0, 10.5, -6.5, 4.5),
)
hidedecorations!(ax)
hidespines!(ax)

# ---------------------------------------------------------
# Geometry: boxes for the two layers
# ---------------------------------------------------------
# Lid: z ∈ (-Hs, 0)  →  y ∈ (0, 3)
# Interior: below     →  y ∈ (-5, 0)
# CMB at y = 3

lid_x = (0.0, 8.0)
lid_y = (0.0, 3.0)
int_y = (-5.0, 0.0)

# Lid region (light blue fill)
poly!(ax, Point2f[(lid_x[1], lid_y[1]), (lid_x[2], lid_y[1]),
                   (lid_x[2], lid_y[2]), (lid_x[1], lid_y[2])];
      color=(:steelblue, 0.12), strokewidth=0)

# Interior region (light orange fill)
poly!(ax, Point2f[(lid_x[1], int_y[1]), (lid_x[2], int_y[1]),
                   (lid_x[2], int_y[2]), (lid_x[1], int_y[2])];
      color=(:coral, 0.10), strokewidth=0)

# Boundaries
lines!(ax, [lid_x[1], lid_x[2]], [3.0, 3.0]; color=:black, linewidth=3.0)     # CMB
lines!(ax, [lid_x[1], lid_x[2]], [0.0, 0.0]; color=:black, linewidth=2.0,
       linestyle=:dash)                                                           # Interface
lines!(ax, [lid_x[1], lid_x[1]], [int_y[1], lid_y[2]]; color=:gray60, linewidth=0.8)
lines!(ax, [lid_x[2], lid_x[2]], [int_y[1], lid_y[2]]; color=:gray60, linewidth=0.8)

# Boundary labels
text!(ax, 4.0, 3.1; text=L"$\mathrm{CMB}\ (z = 0)$", align=(:center, :bottom), fontsize=13)
# text!(ax, 4.3, 0.0; text=L"$\text{Interface}(z = -H_s)$", align=(:left, :center), fontsize=13)

# ---------------------------------------------------------
# Region labels
# ---------------------------------------------------------
# text!(ax, 4.0, 3.7; text=L"$\text{Mantle (rigid)}$",
#       align=(:center, :center), fontsize=13, color=:gray40)

text!(ax, 0.5, 1.5;
      text=L"$\text{Stratified lid} (N = N_s)$",
      align=(:left, :center), fontsize=14, font=:bold, color=:steelblue)

text!(ax, 3.2, -3.0;
      text=L"$\text{Convective interior} (N \approx 0)$",
      align=(:left, :center), fontsize=14, font=:bold, color=(:coral, 0.8))

# ---------------------------------------------------------
# Height annotation: H_s bracket on the left
# ---------------------------------------------------------
bx = -1.0
arrows2d!(ax, [bx], [1.5], [0.0],  [1.5]; color=:black, shaftwidth=1.2, tipwidth=8, tiplength=8)
arrows2d!(ax, [bx], [1.5], [0.0], [-1.5]; color=:black, shaftwidth=1.2, tipwidth=8, tiplength=8)
text!(ax, bx - 0.2, 1.5; text=L"$H_s$", align=(:right, :center), fontsize=15)

# ---------------------------------------------------------
# Background fields: Ω and B₀ arrows on the right
# ---------------------------------------------------------
ax_r = 9.2
# Rotation
arrows2d!(ax, [ax_r], [0.5], [0.0], [2.0]; color=:black, shaftwidth=2.0, tipwidth=12, tiplength=12)
text!(ax, ax_r + 0.15, 1.8; text=L"$\Omega$", align=(:left, :center),
      fontsize=16, color=:black)

# Magnetic field
arrows2d!(ax, [ax_r + 0.8], [0.5], [0.0], [2.0]; color=:purple, shaftwidth=2.0, tipwidth=12, tiplength=12)
text!(ax, ax_r + 0.95, 1.8; text=L"$\mathbf{B}_0$", align=(:left, :center),
      fontsize=16, color=:purple)

# z-axis label
text!(ax, ax_r + 0.4, 2.8; text=L"$\hat{z}$", align=(:center, :bottom), fontsize=14)

# ---------------------------------------------------------
# Incoming fast wave from interior (wiggly arrow)
# ---------------------------------------------------------
x_wave, y_wave = incoming_wave_path()
lines!(ax, x_wave, y_wave; color=:coral, linewidth=2.0)

# Draw only the arrowhead so the incoming-wave marker stays clean at the interface.
incoming_arrowhead = incoming_wave_arrowhead_geometry(x_wave, y_wave)
poly!(ax, Point2f.(incoming_arrowhead.points); color=:coral, strokewidth=0)

text!(ax, 1.3, -3.0; text=L"$\text{Fast wave} (\omega_T)$",
      align=(:center, :center), fontsize=13, color=(:coral, 0.9),
      rotation=π/2)

# ---------------------------------------------------------
# Trapped lid eigenmodes: φ₀, φ₁, φ₂
# ---------------------------------------------------------
# Draw three eigenfunctions cos[(n+1/2)π(z+Hs)/Hs] in the lid
mode_colors = [:navy, :teal, :mediumpurple]
mode_labels = [L"\phi_0", L"\phi_1", L"\phi_2"]
x_offsets   = [5.0, 6.2, 7.2]
amp         = 0.40

nz = 200
zs = range(0.0, 3.0, length=nz)   # mapped lid coordinate

for (ni, (xc, col, lab)) in enumerate(zip(x_offsets, mode_colors, mode_labels))
    n = ni - 1
    # φ_n(z) = cos[(n+½)π(z+Hs)/Hs], mapped to y ∈ [0,3]
    # In plot coords: y=0 is interface (z=-Hs), y=3 is CMB (z=0)
    # φ_n = √2 cos[(n+½)π·(3-y)/3]  (CMB node at y=3)
    ϕ = [sqrt(2) * cos((n + 0.5) * π * y / 3.0) for y in zs]
    xs = xc .+ amp .* ϕ

    lines!(ax, xs, collect(zs); color=col, linewidth=1.8)
    # Zero line
    lines!(ax, [xc, xc], [0.0, 3.0]; color=(col, 0.25), linewidth=0.6, linestyle=:dot)
    # Label at bottom
    text!(ax, xc, -0.10; text=lab, align=(:center, :top), fontsize=14, color=col)
end

text!(ax, 6.1, -0.8; text=L"$\text{Trapped MAC modes}$", align=(:center, :top),
      fontsize=12, color=:gray30)

# ---------------------------------------------------------
# Coupling arrow at interface
# ---------------------------------------------------------
arrows2d!(ax, [3.2], [0.0], [1.3], [0.0]; color=:black, shaftwidth=1.5, tipwidth=10, tiplength=10)
text!(ax, 3.85, 0.25; text=L"$g_n$", align=(:center, :bottom), fontsize=15)

# ---------------------------------------------------------
# Boundary condition annotations (small, near boundaries)
# ---------------------------------------------------------
text!(ax, 0.15, 2.75; text=L"$u_z=0$", align=(:left, :top),
      fontsize=11, color=:gray40)
text!(ax, 0.15, 0.25; text=L"$\partial\phi_n/\partial z=0$",
      align=(:left, :bottom), fontsize=11, color=:gray40)

save("fig1.png", fig, px_per_unit=5)
