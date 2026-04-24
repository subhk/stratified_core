using CairoMakie, Printf

CairoMakie.activate!()

include(joinpath(@__DIR__, "paper_figures_common.jl"))
using .PaperFiguresCommon

set_theme!(Theme(
    fontsize = 20,
    Axis = (
        topspinevisible = false,
        rightspinevisible = false,
        xgridvisible = false,
        ygridvisible = false,
        spinewidth = 1.5,
        xtickwidth = 1.5,
        ytickwidth = 1.5,
        xticklabelsize = 18,
        yticklabelsize = 18,
        xlabelsize = 20,
        ylabelsize = 20,
        xlabelpadding = 6,
        ylabelpadding = 6,
        xtickformat = values -> [@sprintf("%g", v) for v in values],
        ytickformat = values -> [@sprintf("%g", v) for v in values],
    ),
    Lines = (linewidth = 2.0,),
    Legend = (framevisible = false, labelsize = 18),
))

p = ModelParams()
cp = CouplingParams()

omegaT = omegaM(0, p)
times, A, gs, deltas, gammas = rk4_envelope(
    p, cp;
    omegaT = omegaT,
    tmax = 250.0 * YEAR,
    dt = 0.05 * YEAR,
)
ET, En, Etot, epsE = diagnostics(times, A, cp.gammaT, gammas)
AT_exact, B0_exact, g0, gamma0 = exact_single_mode(times, p, cp)

ty = times ./ YEAR
norm0 = ET[1]

fig = Figure(size = (500, 400), fontsize = 20)

ax1 = Axis(
    fig[1, 1],
    xlabel = L"$\text{time (years)}$",
    ylabel = L"$\text{normalized energy}$",
    topspinevisible = false,
    rightspinevisible = false,
)
lines!(ax1, ty, ET ./ norm0; label = L"$E_T / E_T(0)$")
lines!(ax1, ty, En[1] ./ norm0; label = L"$E_0 / E_T(0)$")
lines!(ax1, ty, abs2.(B0_exact) ./ norm0; linestyle = :dash, label = L"$E_0 \text{exact (single-mode)}$")
if length(En) >= 2
    lines!(ax1, ty, En[2] ./ norm0; label = L"$E_1 / E_T(0)$")
end
lines!(ax1, ty, Etot ./ norm0; linestyle = :dashdot, label = L"$E_\text{tot} / E_T(0)$")
axislegend(ax1, position = :rt)

ylims!(ax1, -0.05, 1.05)
xlims!(ax1, -1., 250)
ax1.xticks = 0:50:250
ax1.yticks = 0:0.25:1

# ax2 = Axis(
#     fig[1, 2],
#     xlabel = L"$\text{time (years)}$",
#     ylabel = L"$\text{normalized energy}$",
#     yscale = log10,
#     topspinevisible = false,
#     rightspinevisible = false,
# )
# lines!(ax2, ty, Etot ./ norm0; label = L"$E_\text{tot} / E_T(0)$")
# lines!(ax2, ty, abs.(epsE) .+ 1.0e-20; linestyle = :dash, label = L"$|\epsilon E|$")
# axislegend(ax2, position = :rb)

# ylims!(ax2, 1e-20, 1.05)
# xlims!(ax2, -1, 250)
# ax2.xticks = 0:50:250
# ax2.yticks = [1e-20, 1e-15, 1e-10, 1e-5, 1.0]
# ax2.ytickformat = x -> ["10⁻²⁰", "10⁻¹⁵", "10⁻¹⁰", "10⁻⁵", "1"]

# Label(fig[1, 1, TopLeft()], L"$(a)$", fontsize = 28, padding = (0, 0, 20, 0), halign = :right)
# Label(fig[1, 2, TopLeft()], L"$(b)$", fontsize = 28, padding = (0, 0, 20, 0), halign = :right)

save("fig3.png", fig, px_per_unit=6)

println("Saved fig3.png")
println("Baseline couplings g_n (s^-1) = ", gs)
println("Baseline detunings delta_n (s^-1) = ", deltas)
println("Modal damping rates gamma_n (s^-1) = ", gammas)
