using JLD2, CSVFiles, DataFrames, CairoMakie, ColorSchemes, Colors

function my_theme()
    my_colors = [colorant"#5a0097", colorant"#ac1f71", colorant"#de514c", colorant"#f99c26", colorant"#eefa1b"]
    my_discrete_colors = [colorant"#27f2ee", colorant"#0900ef"]
    colormap = cgrad(ColorScheme(my_colors))
    palette = (color = my_discrete_colors,)
    return Theme(; colormap, palette)
end

##

#region ## Figure Hollow

p = "Output/FigHollow/FigHollow"
files = [
    ["$(p)a.jld2"],
    ["$(p)b.jld2"],
    ["$(p)c.jld2"]]

zscales = [0.7; 0.5; 0.9]

function FigHollow(files, zscales)
    fig = Figure(resolution = (0.9 * 600, 0.9 * 800), font = "CMU Serif Roman")
    local hmap
    for (row, rowfiles) in enumerate(files), (col, file) in enumerate(rowfiles)
        dict = load(file)
        data = sum(values(dict["mjdict"]))
        xs, ys = dict["xs"], real.(dict["ys"])
        m = zscales[row, col] * maximum(data)
        data ./= m
        zlims = (0, 1)

        xlabel = row == 3 ? L"\Phi/\Phi_0" : ""
        ylabel = "ω (meV)"
        label = "LDOS (arb. units)"

        ax = Axis(fig[row, col]; xlabel, ylabel, xlabelpadding = -7)

        row < 3 && hidexdecorations!(ax, ticks = false)
        #col > (row == 1 ? 2 : 1) && hideydecorations!(ax, ticks = false)

        hmap = heatmap!(ax, xs, ys, data; colormap = :thermal, colorrange = zlims)
        cbar = Colorbar(fig, hmap, label = label, ticklabelsvisible = true, ticks = [0,1],
            labelpadding = -10,  width = 10, ticksize = 0,
            ticklabelpad = 5)
        fig[row, col+1] = cbar
    end

    Label(fig[1, 1, TopLeft()], "a", padding = (-40, 0, -20, 0), font = "CMU Serif Bold", textsize = 20)
    Label(fig[2, 1, TopLeft()], "b", padding = (-40, 0, -20, 0), font = "CMU Serif Bold", textsize = 20)
    Label(fig[3, 1, TopLeft()], "c", padding = (-40, 0, -20, 0), font = "CMU Serif Bold", textsize = 20)
    # Label(fig[2, 1, TopLeft()], "d", padding = (-40, 0, -20, 0), font = "CMU Serif Bold", textsize = 20)
    # Label(fig[2, 2, TopLeft()], "e", padding = (35, 0, -20, 0), font = "CMU Serif Bold", textsize = 20)
    # Label(fig[2, 3, TopLeft()], "f", padding = (-20, 0, -20, 0), font = "CMU Serif Bold", textsize = 20)
    # Label(fig[3, 1, TopLeft()], "g", padding = (-40, 0, -20, 0), font = "CMU Serif Bold", textsize = 20)
    # Label(fig[3, 2, TopLeft()], "h", padding = (35, 0, -20, 0), font = "CMU Serif Bold", textsize = 20)
    # Label(fig[3, 3, TopLeft()], "i", padding = (-20, 0, -20, 0), font = "CMU Serif Bold", textsize = 20)
    # Label(fig[1, 1, Top()], L"U_{min}=-30meV", tellwidth=false, tellheight=false, textsize = 25)
    # Label(fig[1, 3, Top()], L"U_{min}=-60meV", tellwidth=false, tellheight=false, textsize = 25)
    # Label(fig[1, 5, Top()], L"U_{min}=-120meV", tellwidth=false, tellheight=false, textsize = 25)

    colgap!(fig.layout, 1, 5)
    rowgap!(fig.layout, 1, 5)
    rowgap!(fig.layout, 2, 5)

    return fig
end

f = with_theme(my_theme()) do
    FigHollow(files, zscales)
end

save("Output/FigHollow/FigHollow_panels.pdf", f)

f

#endregion

#region ## Figure Tubular

p = "Output/FigTubular/FigTubular"
files = [
    ["$(p)b.jld2", "$(p)b.jld2", "$(p)c.jld2"],
    ["$(p)d.jld2", "$(p)e.jld2", "$(p)f.jld2"],
    ["$(p)g.jld2", "$(p)h.jld2", "$(p)i.jld2"]]

zscales = 2 .* [0 1.5 1.3; 1.2 1.2 1.2; 1.1 1.1 1.1]

function FigTubular(files, zscales)
    # fig = Figure(resolution = (1200, 1200), font =:sans)
    fig = Figure(resolution = (1100, 650), font = "CMU Serif Roman")
    local hmap
    for (row, rowfiles) in enumerate(files), (col, file) in enumerate(rowfiles)
        row == col == 1 && continue
        dict = load(file)
        data = sum(values(dict["mjdict"]))
        m = zscales[row, col] * data[end,end]
        data ./= m
        xs, ys = dict["xs"], real.(dict["ys"])
        zlims = (0, 1)

        xlabel = row == 3 ? L"\Phi/\Phi_0" : ""
        ylabel = col > (row == 1 ? 2 : 1) ? "" : "ω (meV)"
        label = "LDOS (arb. units)"

        ax = Axis(fig[row, col]; xlabel, ylabel, xlabelpadding = -7)

        row < 3 && hidexdecorations!(ax, ticks = false)
        col > (row == 1 ? 2 : 1) && hideydecorations!(ax, ticks = false)

        hmap = heatmap!(ax, xs, ys, data; colormap = :thermal, colorrange = zlims)
        if col == 3
            cbar = Colorbar(fig, hmap, label = label, ticklabelsvisible = true, ticks = [0,1],
                labelpadding = -10,  width = 10,ticksize = 2,
                ticklabelpad = 5)
            fig[row, col+1] = cbar
        end
    end

    Label(fig[1, 1, TopLeft()], "a", padding = (-40, 0, -15, 0), font = "CMU Serif Bold", textsize = 20)
    Label(fig[1, 2, TopLeft()], "b", padding = (35, 0, -15, 0), font = "CMU Serif Bold", textsize = 20)
    Label(fig[1, 3, TopLeft()], "c", padding = (-20, 0, -15, 0), font = "CMU Serif Bold", textsize = 20)
    Label(fig[2, 1, TopLeft()], "d", padding = (-40, 0, -20, 0), font = "CMU Serif Bold", textsize = 20)
    Label(fig[2, 2, TopLeft()], "e", padding = (35, 0, -20, 0), font = "CMU Serif Bold", textsize = 20)
    Label(fig[2, 3, TopLeft()], "f", padding = (-20, 0, -20, 0), font = "CMU Serif Bold", textsize = 20)
    Label(fig[3, 1, TopLeft()], "g", padding = (-40, 0, -20, 0), font = "CMU Serif Bold", textsize = 20)
    Label(fig[3, 2, TopLeft()], "h", padding = (35, 0, -20, 0), font = "CMU Serif Bold", textsize = 20)
    Label(fig[3, 3, TopLeft()], "i", padding = (-20, 0, -20, 0), font = "CMU Serif Bold", textsize = 20)
    # Label(fig[1, 1, Top()], L"U_{min}=-30meV", tellwidth=false, tellheight=false, textsize = 25)
    # Label(fig[1, 3, Top()], L"U_{min}=-60meV", tellwidth=false, tellheight=false, textsize = 25)
    # Label(fig[1, 5, Top()], L"U_{min}=-120meV", tellwidth=false, tellheight=false, textsize = 25)

    colgap!(fig.layout, 1, -30)
    colgap!(fig.layout, 2, 25)
    colgap!(fig.layout, 3, 5)
    rowgap!(fig.layout, 1, 5)
    rowgap!(fig.layout, 2, 5)

    return fig
end

f = with_theme(my_theme()) do
    FigTubular(files, zscales)
end

save("Output/FigTubular/FigTubular_panels.pdf", f)

f

#endregion

#region ## Figure WF

p = "Output/FigWF/FigWF"
files = [
    ["$(p)a.jld2", "$(p)b.jld2"],
    ["$(p)c.jld2", "$(p)d.jld2"]]

zscales = [0.8 0.8; 0.7 0.7]

function FigWF(files, zscales)
    fig = Figure(resolution = (700, 500), font = "CMU Serif Roman")
    local hmap
    for (row, rowfiles) in enumerate(files), (col, file) in enumerate(rowfiles)
        dict = load(file)
        data = sum(values(dict["mjdict"]))
        xs, ys = dict["xs"], real.(dict["ys"])
        m = zscales[row, col] * maximum(data)
        data ./= m
        zlims = (0, 1)

        xlabel = row == 2 ? L"\Phi/\Phi_0" : ""
        ylabel = "ω (meV)"
        label = "LDOS (arb. units)"

        ax = Axis(fig[row, col]; xlabel, ylabel, xlabelpadding = -7)

        row < 2 && hidexdecorations!(ax, ticks = false)
        col > 1 && hideydecorations!(ax, ticks = false)

        hmap = heatmap!(ax, xs, ys, data; colormap = :thermal, colorrange = zlims)
        if col == 2
            cbar = Colorbar(fig, hmap, label = label, ticklabelsvisible = true, ticks = [0,1],
                labelpadding = -10,  width = 10,ticksize = 2,
                ticklabelpad = 5)
            fig[row, col+1] = cbar
        end
    end

    Label(fig[1, 1, TopLeft()], "a", padding = (-40, 0, -20, 0), font = "CMU Serif Bold", textsize = 20)
    Label(fig[1, 2, TopLeft()], "b", padding = (-20, 0, -20, 0), font = "CMU Serif Bold", textsize = 20)
    Label(fig[2, 1, TopLeft()], "c", padding = (-40, 0, -20, 0), font = "CMU Serif Bold", textsize = 20)
    Label(fig[2, 2, TopLeft()], "d", padding = (-20, 0, -20, 0), font = "CMU Serif Bold", textsize = 20)
    # Label(fig[2, 2, TopLeft()], "e", padding = (35, 0, -20, 0), font = "CMU Serif Bold", textsize = 20)
    # Label(fig[2, 3, TopLeft()], "f", padding = (-20, 0, -20, 0), font = "CMU Serif Bold", textsize = 20)
    # Label(fig[3, 1, TopLeft()], "g", padding = (-40, 0, -20, 0), font = "CMU Serif Bold", textsize = 20)
    # Label(fig[3, 2, TopLeft()], "h", padding = (35, 0, -20, 0), font = "CMU Serif Bold", textsize = 20)
    # Label(fig[3, 3, TopLeft()], "i", padding = (-20, 0, -20, 0), font = "CMU Serif Bold", textsize = 20)
    # Label(fig[1, 1, Top()], L"U_{min}=-30meV", tellwidth=false, tellheight=false, textsize = 25)
    # Label(fig[1, 3, Top()], L"U_{min}=-60meV", tellwidth=false, tellheight=false, textsize = 25)
    # Label(fig[1, 5, Top()], L"U_{min}=-120meV", tellwidth=false, tellheight=false, textsize = 25)

    colgap!(fig.layout, 1, 25)
    colgap!(fig.layout, 2, 5)
    rowgap!(fig.layout, 1, 5)

    return fig
end

f = with_theme(my_theme()) do
    FigWF(files, zscales)
end

save("Output/FigWF/FigWF_panels.pdf", f)

f

#endregion

#region ## Figure Solid

p = "Output/FigSolid/FigSolid_"
files = [
    ["$(p)A_ldos.jld2", "$(p)A_cond_short.jld2", "$(p)A_cond_long.jld2"],
    ["$(p)B_ldos.jld2", "$(p)B_cond_short.jld2", "$(p)B_cond_long.jld2"],
    ["$(p)C_ldos.jld2", "$(p)C_cond_short.jld2", "$(p)C_cond_long.jld2"]]

# zscales = [0.5 0.4 0.4; 0.44 0.4 0.4; 0.3 0.25 0.4]
zscales = 2 .* [0.8 0.8 0.8; 0.9 0.9 0.9; 1 1 1]

function FigSolid(files, zscales)
    # fig = Figure(resolution = (1200, 1200), font =:sans)
    fig = Figure(resolution = (1100, 650), font = "CMU Serif Roman")
    local hmap
    for (col, colfiles) in enumerate(files), (row, file) in enumerate(colfiles)
        dict = load(file)
        data = sum(values(dict["mjdict"]))
        # maxref = zscales[row, col] * maximum(data)
        maxref = zscales[row, col] * data[end,end]
        m = row == 1 ? maxref : 1.0
        data ./= m
        xs, ys = dict["xs"], real.(dict["ys"])
        zlims = row == 1 ? (0, 1) : (0, maxref)

        xlabel = row == 3 ? L"\Phi/\Phi_0" : ""
        ylabel = col == 1 ? (row == 1 ? "ω (meV)" : "V (mV)") : ""
        label = col == 3 ? row == 1 ? "LDOS (arb. units)" : L"dI/dV\,\,(G_0)" : ""

        ax = Axis(fig[row, 2*col-1]; xlabel, ylabel, xlabelpadding = -7)

        row < 3 && hidexdecorations!(ax, ticks = false)
        col > 1 && hideydecorations!(ax, ticks = false)

        hmap = heatmap!(ax, xs, ys, data; colormap = :thermal, colorrange = zlims)
        cbar = Colorbar(fig, hmap, label = label, ticklabelsvisible = true, ticks = ifelse(row == 1, [0,1], Makie.automatic),
            labelpadding = 5,  width = 10,ticksize = 2,
            ticklabelpad = 5)
        fig[row, 2*col] = cbar
    end

    Label(fig[1, 1, TopLeft()], "b", padding = (-40, 0, 10, 0), font = "CMU Serif Bold", textsize = 20)
    Label(fig[1, 3, TopLeft()], "c", padding = (-20, 0, 10, 0), font = "CMU Serif Bold", textsize = 20)
    Label(fig[1, 5, TopLeft()], "d", padding = (-20, 0, 10, 0), font = "CMU Serif Bold", textsize = 20)
    Label(fig[2, 1, TopLeft()], "f", padding = (-40, 0, -20, 0), font = "CMU Serif Bold", textsize = 20)
    Label(fig[2, 3, TopLeft()], "g", padding = (-20, 0, -20, 0), font = "CMU Serif Bold", textsize = 20)
    Label(fig[2, 5, TopLeft()], "h", padding = (-20, 0, -20, 0), font = "CMU Serif Bold", textsize = 20)
    Label(fig[3, 1, TopLeft()], "j", padding = (-40, 0, -20, 0), font = "CMU Serif Bold", textsize = 20)
    Label(fig[3, 3, TopLeft()], "k", padding = (-20, 0, -20, 0), font = "CMU Serif Bold", textsize = 20)
    Label(fig[3, 5, TopLeft()], "l", padding = (-20, 0, -20, 0), font = "CMU Serif Bold", textsize = 20)
    Label(fig[1, 1, Top()], L"U_{min}=-40meV", tellwidth=false, tellheight=false, textsize = 20)
    Label(fig[1, 3, Top()], L"U_{min}=-70meV", tellwidth=false, tellheight=false, textsize = 20)
    Label(fig[1, 5, Top()], L"U_{min}=-140meV", tellwidth=false, tellheight=false, textsize = 20)

    colgap!(fig.layout, 1, 5)
    colgap!(fig.layout, 3, 5)
    colgap!(fig.layout, 5, 5)
    rowgap!(fig.layout, 1, 5)
    rowgap!(fig.layout, 2, 5)

    return fig
end

f = with_theme(my_theme()) do
    FigSolid(files, zscales)
end

save("Output/FigSolid/FigSolid_panels.pdf", f)

f

#endregion

#region ## Figure Effective

p = "Output/FigEffective/FigEffective"
files = [["$(p)_a_ldos.jld2", "$(p)_b_ldos.jld2"],
         ["$(p)_c_ldos.jld2", "$(p)_d_ldos.jld2"]]

zscales = [0.6 0.4; 0.6 0.4]

function FigEffective(files, zscales)
    # fig = Figure(resolution = (1200, 1200), font =:sans)
    fig = Figure(resolution = (0.9 * 600, 0.9 * 800), font = "CMU Serif Roman")
    local hmap
    for (row, rowfiles) in enumerate(files), (col, file) in enumerate(rowfiles)
        dict = load(file)

        data = sum(values(dict["mjdict"]))
        m = zscales[row, col] * maximum(data)
        data ./= m
        xs, ys = dict["xs"], real.(dict["ys"])
        zlims = (0, 1)

        if col == 2
            data ./= xs
            data .*= maximum(xs)
            data .= sqrt.(data)
        end
        data = [reverse(data[:, 2:end], dims = 2) data]
        ys = [-reverse(ys[2:end]); ys]

        xlabel = col == 1 ? L"\Phi/\Phi_0" : "r (nm)"
        ylabel = col == 1 ? "ω (meV)" : ""
        label =  "LDOS (arb. units)"
        xticks = col == 2 ? (0:20:60) : (0:3)

        ax = Axis(fig[row, col]; xlabel, ylabel, xticks)

        # row < 3 && hidexdecorations!(ax, ticks = false)
        col > 1 && hideydecorations!(ax, ticks = false)

        hmap = heatmap!(ax, xs, ys, data; colormap = :thermal, colorrange = zlims)
        row == 1 && hidexdecorations!(ax, ticks = false)
        if col == 1
            text!(ax,
                (2.4, 0.23),
                text = (row == 1 ? "Solid core" : "Modified hollow core"),
                color = (:white, 0.99),
                align = (:right, :center))
            if row == 2
                for (pos, text) in zip(
                        [(1.25,0.185), (1.25,0.15), (1.25,0.11), (1.25,0.07), (1.25,0.025), (1.25,-0.185), (1.25, -0.15), (1.25, -0.11), (1.25, -0.07), (1.25, -0.025)], 
                        ["-1/2", "1/2", "3/2", "5/2", "7/2","1/2", "-1/2", "-3/2", "-5/2", "-7/2"])
                    text!(ax, pos; text, color = (:white, 0.99), align = (:right, :center), textsize = 15)
                end
            end
        else
            xlims!(ax, (0, 70))
            ax.backgroundcolor = :black
            cbar = Colorbar(fig, hmap, label = label, ticklabelsvisible = true, ticks = [0,1],
                labelpadding = -10,  width = 10, ticksize = 2,
                ticklabelpad = 5)
            fig[row, 3] = cbar
        end

        # Dashed vertical line
        if (row, col) != (2, 2)
            linex = col == 1 ? 1 : 56
            linecolor = col == 1 ? RGBA(1,1,1,0.9) : RGBA(0,0,0,0.9)
            lines!(ax, [linex,linex], [-0.26, 0.26], color = linecolor, linewidth = 1, linestyle = [10,14])
        end

        # Rav on top
        if (row, col) == (1, 2)
            Rav = 56
            ax2 = Axis(fig[row, col], xaxisposition = :top, xticks = ([Rav], [L"R_\textrm{av}"]))
            ax2.backgroundcolor = :black
            xlims!(ax2, extrema(xs))
            hidespines!(ax2)
            hideydecorations!(ax2)
        end
    end

    Label(fig[1, 1, TopLeft()], "a", padding = (-40, 0, -20, 0), font = "CMU Serif Bold", textsize = 20)
    Label(fig[1, 2, TopLeft()], "b", padding = (-20, 0, -20, 0), font = "CMU Serif Bold", textsize = 20)
    Label(fig[2, 1, TopLeft()], "c", padding = (-40, 0, -20, 0), font = "CMU Serif Bold", textsize = 20)
    Label(fig[2, 2, TopLeft()], "d", padding = (-20, 0, -20, 0), font = "CMU Serif Bold", textsize = 20)

    colsize!(fig.layout, 2, Aspect(1, 0.25))
    colgap!(fig.layout, 1, 23)
    colgap!(fig.layout, 2, 10)

    return fig
end

f = with_theme(my_theme()) do
    FigEffective(files, zscales)
end

save("Output/FigEffective/FigEffective.pdf", f)

f

#endregion

#region ## Figure SOC

p = "Output/FigSOC"
files = [
    ["$p/FigNOSOC_notop_ldos.jld2", "$p/FigSOC_notop_ldos.jld2"],
    ["$p/FigNOSOC_top_ldos.jld2", "$p/FigSOC_top_ldos.jld2"]]

zscaledos = 2

function FigSOC(files)
    # fig = Figure(resolution = (1200, 1200), font =:sans)
    fig = Figure(resolution = (600, 500), font = "CMU Serif Roman")
    local hmap
    for (row, rowfiles) in enumerate(files), (col, file) in enumerate(rowfiles)
        dict = load(file)

        data = sum(values(dict["mjdict"]))
        m = zscaledos * data[end, end]
        data ./= m
        xs, ys = dict["xs"], real.(dict["ys"])
        zlims = (0, 1)

        data = [reverse(data[:, 2:end], dims = 2) data]
        ys = [-reverse(ys[2:end]); ys]

        xlabel = row == 2 ? L"\Phi/\Phi_0" : ""
        ylabel = col == 1 ? "ω (meV)" : ""
        labelbar =  "LDOS (arb. units)"

        ax = Axis(fig[row, col]; xlabel, ylabel)

        row < 2 && hidexdecorations!(ax, ticks = false)
        col > 1 && hideydecorations!(ax, ticks = false)

        hmap = heatmap!(ax, xs, ys, data; colormap = :thermal, colorrange = zlims, rasterize = 5)
        # if col == 2
        #     cbar = Colorbar(fig, hmap, label = labelbar, ticklabelsvisible = false, labelpadding = 5, width = 10, ticksize = 2,
        #     ticklabelpad = 5)
        #     fig[row, 3] = cbar
        # end
        if (row, col) == (2, 2)
            text!(ax, (1.15, 0.00); text = "Majorana", color = (:white, 0.99), align = (:center, :bottom), textsize = 15)
        end
    end

    Label(fig[1, 1, TopLeft()], "a", padding = (-40, 0, 0, 0), font = "CMU Serif Bold", textsize = 20)
    Label(fig[1, 2, TopLeft()], "c", padding = (-20, 0, 0, 0), font = "CMU Serif Bold", textsize = 20)
    Label(fig[2, 1, TopLeft()], "b", padding = (-40, 0, 0, 0), font = "CMU Serif Bold", textsize = 20)
    Label(fig[2, 2, TopLeft()], "d", padding = (-20, 0, 0, 0), font = "CMU Serif Bold", textsize = 20)
    titles = ["α = g = 0" "α, g ≠ 0, non-topological";
             "α = g = 0" "α, g ≠ 0, topological"]
    for row in 1:2, col in 1:2
        Label(fig[row, col, Top()], titles[row, col], font = "CMU Serif Italic", textsize = 18, tellwidth=false, tellheight=true)
    end

    colgap!(fig.layout, 1, 20)
    rowgap!(fig.layout, 1, 5)

    return fig
end

f = with_theme(my_theme()) do
    FigSOC(files)
end

save("Output/FigSOC/FigSOC.pdf", f)

f

#endregion

#region ## Figure Tubular Appendix

p = "Output/FigTubularAppendix/FigTubular"
files = [
    ["$(p)b.jld2", "$(p)b.jld2", "$(p)c.jld2"],
    ["$(p)d.jld2", "$(p)e.jld2", "$(p)f.jld2"],
    ["$(p)g.jld2", "$(p)h.jld2", "$(p)i.jld2"]]

zscales = 2 .* [0 1.5 1.3; 1.2 1.2 1.2; 1.1 1.1 1.1]

function FigTubular(files, zscales)
    # fig = Figure(resolution = (1200, 1200), font =:sans)
    fig = Figure(resolution = (1100, 650), font = "CMU Serif Roman")
    # fig = Figure(resolution = (1100, 650), font = "CMU Serif Roman")
    local hmap
    for (row, rowfiles) in enumerate(files), (col, file) in enumerate(rowfiles)
        row == col == 1 && continue
        dict = load(file)
        data = sum(values(dict["mjdict"]))
        m = zscales[row, col] * data[end,end]
        data ./= m
        xs, ys = dict["xs"], real.(dict["ys"])
        zlims = (0, 1)

        xlabel = row == 3 ? L"\Phi/\Phi_0" : ""
        ylabel = col > (row == 1 ? 2 : 1) ? "" : "ω (meV)"
        label = "LDOS (arb. units)"

        ax = Axis(fig[row, col]; xlabel, ylabel, xlabelpadding = -7)

        row < 3 && hidexdecorations!(ax, ticks = false)
        col > (row == 1 ? 2 : 1) && hideydecorations!(ax, ticks = false)

        hmap = heatmap!(ax, xs, ys, data; colormap = :thermal, colorrange = zlims)
        if col == 3
            cbar = Colorbar(fig, hmap, label = label, ticklabelsvisible = true, ticks = [0,1],
                labelpadding = -10,  width = 10,ticksize = 2,
                ticklabelpad = 5)
            fig[row, col+1] = cbar
        end
        for linex in [0.5, 1.5, 2.5]
            linecolor = RGBA(1,1,1,0.9)
            lines!(ax, [linex,linex], [-0.26, 0.26], color = linecolor, linewidth = 1, linestyle = [10,14])
        end
    end

    Label(fig[1, 1, TopLeft()], "a", padding = (-40, 0, -15, 0), font = "CMU Serif Bold", textsize = 20)
    Label(fig[1, 2, TopLeft()], "b", padding = (35, 0, -15, 0), font = "CMU Serif Bold", textsize = 20)
    Label(fig[1, 3, TopLeft()], "c", padding = (-20, 0, -15, 0), font = "CMU Serif Bold", textsize = 20)
    Label(fig[2, 1, TopLeft()], "d", padding = (-40, 0, -20, 0), font = "CMU Serif Bold", textsize = 20)
    Label(fig[2, 2, TopLeft()], "e", padding = (35, 0, -20, 0), font = "CMU Serif Bold", textsize = 20)
    Label(fig[2, 3, TopLeft()], "f", padding = (-20, 0, -20, 0), font = "CMU Serif Bold", textsize = 20)
    Label(fig[3, 1, TopLeft()], "g", padding = (-40, 0, -20, 0), font = "CMU Serif Bold", textsize = 20)
    Label(fig[3, 2, TopLeft()], "h", padding = (35, 0, -20, 0), font = "CMU Serif Bold", textsize = 20)
    Label(fig[3, 3, TopLeft()], "i", padding = (-20, 0, -20, 0), font = "CMU Serif Bold", textsize = 20)
    # Label(fig[1, 1, Top()], L"U_{min}=-30meV", tellwidth=false, tellheight=false, textsize = 25)
    # Label(fig[1, 3, Top()], L"U_{min}=-60meV", tellwidth=false, tellheight=false, textsize = 25)
    # Label(fig[1, 5, Top()], L"U_{min}=-120meV", tellwidth=false, tellheight=false, textsize = 25)

    colgap!(fig.layout, 1, -30)
    colgap!(fig.layout, 2, 25)
    colgap!(fig.layout, 3, 5)
    rowgap!(fig.layout, 1, 5)
    rowgap!(fig.layout, 2, 5)

    return fig
end

f = with_theme(my_theme()) do
    FigTubular(files, zscales)
end

save("Output/FigTubularAppendix/FigTubular_panels.pdf", f)

f

#endregion

#region ## Figure Solid Appendix

p = "Output/FigSolidAppendix/FigSolid_"
files = [
    ["$(p)A_ldos.jld2", "$(p)A_cond_short.jld2", "$(p)A_cond_long.jld2"],
    ["$(p)B_ldos.jld2", "$(p)B_cond_short.jld2", "$(p)B_cond_long.jld2"],
    ["$(p)C_ldos.jld2", "$(p)C_cond_short.jld2", "$(p)C_cond_long.jld2"]]

# zscales = [0.5 0.4 0.4; 0.44 0.4 0.4; 0.3 0.25 0.4]
zscales = 2 .* [0.8 0.8 0.8; 0.9 0.9 0.9; 1 1 1]

function FigSolid(files, zscales)
    # fig = Figure(resolution = (1200, 1200), font =:sans)
    fig = Figure(resolution = (1100, 650), font = "CMU Serif Roman")
    local hmap
    for (col, colfiles) in enumerate(files), (row, file) in enumerate(colfiles)
        dict = load(file)
        data = sum(values(dict["mjdict"]))
        # maxref = zscales[row, col] * maximum(data)
        maxref = zscales[row, col] * data[end,end]
        m = row == 1 ? maxref : 1.0
        data ./= m
        xs, ys = dict["xs"], real.(dict["ys"])
        zlims = row == 1 ? (0, 1) : (0, maxref)

        xlabel = row == 3 ? L"\Phi/\Phi_0" : ""
        ylabel = col == 1 ? (row == 1 ? "ω (meV)" : "V (mV)") : ""
        label = col == 3 ? row == 1 ? "LDOS (arb. units)" : L"dI/dV\,\,(G_0)" : ""

        ax = Axis(fig[row, 2*col-1]; xlabel, ylabel, xlabelpadding = -7)

        row < 3 && hidexdecorations!(ax, ticks = false)
        col > 1 && hideydecorations!(ax, ticks = false)

        hmap = heatmap!(ax, xs, ys, data; colormap = :thermal, colorrange = zlims)
        cbar = Colorbar(fig, hmap, label = label, ticklabelsvisible = true, ticks = ifelse(row == 1, [0,1], Makie.automatic),
            labelpadding = 5,  width = 10,ticksize = 2,
            ticklabelpad = 5)
        fig[row, 2*col] = cbar
    end

    Label(fig[1, 1, TopLeft()], "b", padding = (-40, 0, 10, 0), font = "CMU Serif Bold", textsize = 20)
    Label(fig[1, 3, TopLeft()], "c", padding = (-20, 0, 10, 0), font = "CMU Serif Bold", textsize = 20)
    Label(fig[1, 5, TopLeft()], "d", padding = (-20, 0, 10, 0), font = "CMU Serif Bold", textsize = 20)
    Label(fig[2, 1, TopLeft()], "f", padding = (-40, 0, -20, 0), font = "CMU Serif Bold", textsize = 20)
    Label(fig[2, 3, TopLeft()], "g", padding = (-20, 0, -20, 0), font = "CMU Serif Bold", textsize = 20)
    Label(fig[2, 5, TopLeft()], "h", padding = (-20, 0, -20, 0), font = "CMU Serif Bold", textsize = 20)
    Label(fig[3, 1, TopLeft()], "j", padding = (-40, 0, -20, 0), font = "CMU Serif Bold", textsize = 20)
    Label(fig[3, 3, TopLeft()], "k", padding = (-20, 0, -20, 0), font = "CMU Serif Bold", textsize = 20)
    Label(fig[3, 5, TopLeft()], "l", padding = (-20, 0, -20, 0), font = "CMU Serif Bold", textsize = 20)
    Label(fig[1, 1, Top()], L"U_{min}=-40meV", tellwidth=false, tellheight=false, textsize = 20)
    Label(fig[1, 3, Top()], L"U_{min}=-70meV", tellwidth=false, tellheight=false, textsize = 20)
    Label(fig[1, 5, Top()], L"U_{min}=-140meV", tellwidth=false, tellheight=false, textsize = 20)

    colgap!(fig.layout, 1, 5)
    colgap!(fig.layout, 3, 5)
    colgap!(fig.layout, 5, 5)
    rowgap!(fig.layout, 1, 5)
    rowgap!(fig.layout, 2, 5)

    return fig
end

f = with_theme(my_theme()) do
    FigSolid(files, zscales)
end

save("Output/FigSolidAppendix/FigSolid_panels.pdf", f)

f

#endregion

#region ## Figure SOC Appendix

p = "Output/FigSOCAppendix"
files = [
    ["$p/FigNOSOC_notop_ldos.jld2", "$p/FigSOC_notop_ldos.jld2"]]

zscaledos = 2

function FigSOC(files)
    # fig = Figure(resolution = (1200, 1200), font =:sans)
    fig = Figure(resolution = (600, 250), font = "CMU Serif Roman")
    local hmap
    for (row, rowfiles) in enumerate(files), (col, file) in enumerate(rowfiles)
        dict = load(file)

        data = sum(values(dict["mjdict"]))
        m = zscaledos * data[end, end]
        data ./= m
        xs, ys = dict["xs"], real.(dict["ys"])
        zlims = (0, 1)

        data = [reverse(data[:, 2:end], dims = 2) data]
        ys = [-reverse(ys[2:end]); ys]

        xlabel = L"\Phi/\Phi_0"
        ylabel = col == 1 ? "ω (meV)" : ""
        labelbar =  "LDOS (arb. units)"

        ax = Axis(fig[row, col]; xlabel, ylabel)

        # row < 2 && hidexdecorations!(ax, ticks = false)
        col > 1 && hideydecorations!(ax, ticks = false)

        hmap = heatmap!(ax, xs, ys, data; colormap = :thermal, colorrange = zlims, rasterize = 5)
        if col == 2
            cbar = Colorbar(fig, hmap, label = labelbar, ticklabelsvisible = false, labelpadding = 5, width = 10, ticksize = 2,
            ticklabelpad = 5)
            fig[row, 3] = cbar
        end
        # if (row, col) == (2, 2)
        #     text!(ax, (1.15, 0.00); text = "Majorana", color = (:white, 0.99), align = (:center, :bottom), textsize = 15)
        # end
    end

    Label(fig[1, 1, TopLeft()], "a", padding = (-40, 0, 0, 0), font = "CMU Serif Bold", textsize = 20)
    Label(fig[1, 2, TopLeft()], "b", padding = (-20, 0, 0, 0), font = "CMU Serif Bold", textsize = 20)
    titles = ["α = g = 0" "α, g ≠ 0"]
    for row in 1:1, col in 1:2
        Label(fig[row, col, Top()], titles[row, col], font = "CMU Serif Italic", textsize = 18, tellwidth=false, tellheight=true)
    end

    colgap!(fig.layout, 1, 20)
    colgap!(fig.layout, 2, 5)
    # rowgap!(fig.layout, 1, 5)

    return fig
end

f = with_theme(my_theme()) do
    FigSOC(files)
end

save("Output/FigSOCAppendix/FigSOC.pdf", f)

f

#endregion
