using CairoMakie, Printf, LaTeXStrings

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

hs_km_list = [50.0, 100.0, 150.0]
hs_list = hs_km_list .* 1.0e3
modes = collect(0:8)
modes_res = collect(0:8)

fig = Figure(size = (1000, 400), fontsize=34)

ax1 = Axis(
    fig[1, 1],
    xlabel = L"\text{mode} $n$",
    ylabel = L"\text{period (years)}",
    topspinevisible = false,
    rightspinevisible = false,
)

for (Hs_km, Hs) in zip(hs_km_list, hs_list)
    p = ModelParams(Hs = Hs)
    periods = [2.0 * pi / omegaM(n, p) / YEAR for n in modes]
    hs_int = Int(round(Hs_km))
    scatterlines!(ax1, modes, periods; marker = :circle, label = latexstring("H_s = $hs_int \\, \\mathrm{km}"))
end
axislegend(ax1, position = :lb)

ylims!(ax1, 0.0, 13)
xlims!(ax1, -0.1, 8.1)
ax1.xticks = 0:2:8
ax1.yticks = 0:3:12

ax2 = Axis(
    fig[1, 2],
    xlabel = L"\text{mode} $n$",
    ylabel = L"$(\omega_{M,n+1} - \omega_{M,n}) / \gamma_n$",
    topspinevisible = false,
    rightspinevisible = false,
)

for (Hs_km, Hs) in zip(hs_km_list, hs_list)
    p = ModelParams(Hs = Hs)
    ratios = [mode_spacing_ratio(n, p) for n in modes_res]
    hs_int = Int(round(Hs_km))
    scatterlines!(ax2, modes_res, ratios; marker = :rect,
            label = latexstring("H_s = $hs_int \\, \\mathrm{km}"))
end
#hlines!(ax2, [1.0]; linestyle = :dash, color = :black, linewidth = 1.5, label = L"spacing $=$ damping")

axislegend(ax2, position = :rc)

xlims!(ax2, -0.1, 8.1)
ylims!(ax2,  0.0, 1.21)
ax2.xticks = 0:2:8
ax2.yticks = 0:0.4:1.2

Label(fig[1, 1, TopLeft()], L"$(a)$", fontsize = 28, padding = (0, 0, 10, 0), halign = :right)
Label(fig[1, 2, TopLeft()], L"$(b)$", fontsize = 28, padding = (0, 0, 10, 0), halign = :right)

save("fig2.png", fig, pt_per_unit = 6)

