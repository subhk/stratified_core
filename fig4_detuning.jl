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
    ScatterLines = (linewidth = 2.0, markersize = 9),
    Legend = (framevisible = false, labelsize = 18),
))

p = ModelParams()
cp = CouplingParams()

deltas = collect(range(-4.0e-9, 4.0e-9, length = 81))
C_exact = zeros(length(deltas))
C_multi = zeros(length(deltas))

for j in eachindex(deltas)
    delta = deltas[j]
    C_exact[j] = peak_converted_fraction_exact(
        p, cp;
        delta = delta,
        tmax = 250.0 * YEAR,
        dt = 0.10 * YEAR,
    )
    C_multi[j] = peak_converted_fraction_multimode(
        p, cp;
        delta = delta,
        tmax = 250.0 * YEAR,
        dt = 0.10 * YEAR,
    )
end

x = deltas .* YEAR  # rad / yr, since delta has units s^-1

fig = Figure(size = (600, 400), fontsize = 18)
ax = Axis(
    fig[1, 1],
    xlabel = L"$\text{detuning } \delta \text{ (rad yr}^{-1}\text{)}$",
    ylabel = L"$\mathcal{R}_{0,\mathrm{peak}}$",
    topspinevisible = false,
    rightspinevisible = false,
)

scatterlines!(ax, x, C_multi; marker = :circle, label = L"$\text{RK4 multimode}$")
lines!(ax, x, C_exact; linestyle = :dash, label = L"$\text{exact single-mode}$")
vlines!(ax, [0.0]; linestyle = :dot, color = :black, linewidth = 1.5,)
axislegend(ax, position = :rt)

#Label(fig[1, 1, TopLeft()], L"$(a)$", fontsize = 28, padding = (0, 0, 20, 0), halign = :right)

save("fig4.png", fig, px_per_unit=6)

println("Saved fig4.png")
