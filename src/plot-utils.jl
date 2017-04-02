import PLplot
import KernelDensity # for violin plot

# function distinguishable(len, c=colorant"white")
#     cm = Colors.distinguishable_colors(len,c)
#     return hcat([ Cint[c.r.i, c.g.i, c.b.i] for c in cm]...)'
# end

# function distinguishable2(len, c=colorant"white")
#     cm = Colors.distinguishable_colors(len, [LCHab(70, 60, 240)],
#                                     transform=c -> deuteranopic(c, 0.5),
#                                     lchoices=Float64[65, 70, 75, 80],
#                                     cchoices=Float64[0, 50, 60, 70],
#                                     hchoices=linspace(0, 330, 24))
#     cmrgb = map(c->convert(RGB{U8}, c), cm)
#     return hcat([ Cint[c.r.i, c.g.i, c.b.i] for c in cmrgb]...)'
# end

# function colormap(len, cn="blues")
#     cm = convert(Array{RGB{U8},1}, Colors.colormap(cn, len))
#     return hcat([ Cint[c.r.i, c.g.i, c.b.i] for c in cm]...)'
# end

# """Gadfly-like colors: 0-white, 1-black, 2...-Gadfly colors"""
# function setcolormap(len)
#     rgbcm = vcat(fill(Cint(255),3)', zeros(Cint,3)', distinguishable2(len+1))
#     PLplot.setcolormap(rgbcm, index=0)
# end

function fillcoord(lat, lon, col)
    PLplot.plcol1(col)
    xs = [lon-180, lon-180, lon+.999-180, lon+.999-180]
    ys = [lat-90, lat+.999-90, lat+.999-90, lat-90]
    PLplot.plfill(length(xs), xs, ys)
end

function clplot(coords; dev="xwin", fname="", xres=0, yres=0)
    miny = -65
    maxy = 90
    minx = -180
    maxx = minx + 360

    # Setup colors
    lcl = map(c->c[3], coords) |> maximum
    # rgbcm = distinguishable(lcl,colorant"white")
    # PLplot.setcolormap(rgbcm, index=1)
    # PLplot.plscolbg(255,255,255)
    # PLplot.plscol0(1, 0, 0, 0)

    PLplot.plsdev(dev)
    !isempty(fname) && plsfnam(fname) # set file name
    xres != 0 && yres != 0 && PLplot.pageparams!(xlen=xres, ylen=yres) # set resolution
    PLplot.plinit()
    PLplot.pladv(0)
    PLplot.plvpor( .05, .98, .08, .95 )
    PLplot.plwind( minx, maxx, miny, maxy )
    PLplot.plbox( "bcnst", 0.0, 0, "bcnstv", 0.0, 0 )
    # plenv( minx, maxx, miny, maxy, Cint(PLplot.Independent), Cint(PLplot.Default) )

    for (lat, long, col) in coords
        # println(lat, ", ", long, " , ", col)
        fillcoord(lat, long, col/lcl)
    end

    PLplot.plend()
end

# Spread of the clusters in  clustering

function spread(cls, X, I; pdev=:xwin, rows = 4, cols = 4, filename="")
    # pcols = Cint(5)
    # prows = ceil(Cint,(length(cls)-1)/pcols)
    # pcols, prows = 4,4

    pdev, fn = if isempty(filename)
        :xwin, ""
    else
        :pdfcairo, filename
    end

    PLplot.plscolbg(255,255,255)
    PLplot.plscol0(1, 0, 0, 0)
    PLplot.draw(pdev, filename=fn) do opts
        for (i,m) in enumerate(cls[1:end-1])
            i%(cols*rows) == 1 && PLplot.plssub(cols, rows)
            PLplot.histogram(LMCLUS.project(m, clusterdata(X, I, m))[1,:])
            PLplot.labels("spread", "", "Cluster $i")
            if i%(cols*rows) == 0
                PLplot.pleop()
                PLplot.plbop()
            end
        end
    end
end

# Parallel coordinates plot

function parcoord{T<:Real}(X::Matrix{T})

end


# Violin plot

function violin_coords(ys, trim, npoints, bw)
    window = npoints > 1 ? (bw == 0.0 ? KernelDensity.default_bandwidth(ys) : bw) : 0.1
    kd = KernelDensity.kde(ys, bandwidth=window, npoints=npoints)
    idxs = if trim
        ymin, ymax = extrema(ys)
        inside = Bool[ ymin <= d <= ymax for d in kd.x]
        return (kd.density[inside], collect(Float64, kd.x[inside]))
    end
    kd.density, collect(Float64, kd.x)
end

function violin(lbls::Vector, data::Vector{Vector{PLplot.PLFLT}}; trim=false, npoints=200, bw=0.0, width=0.8)
    @assert length(lbls) == length(data) "Dimension of labels and data arrays should be of the same size"

    dim = length(data)
    xlbls = convert(Vector{PLplot.PLFLT}, lbls)

    # Set initial plot dimensions
    ymin = Inf
    ymax = 0.
    xmin, xmax = extrema(xlbls)
    xw = (xmax-xmin)/(length(xlbls)-1)
    xmin -= xw
    xmax += xw

    # Calculate violin shapes per column of input matrix
    vlnw = fill(Float64[],dim)
    vlnc = fill(Float64[],dim)
    for i in 1:dim
        xc = Float64(i)
        widths, centers = violin_coords(data[i], trim, npoints, bw)
        if length(widths) > 0
            tmin, tmax = extrema(centers)
            ymin = min(ymin, tmin)
            ymax = max(ymax, tmax)
            vlnc[i] = centers

            maxwidth = maximum(widths)
            broadcast!(*, widths, widths, 0.5 * width)
            broadcast!(/, widths, widths, tmax)
            vlnw[i] = widths
        else
            vlnw[i] = [median(data[i])]
        end
    end

    PLplot.plenv(xmin, xmax, ymin-1.3, ymax, Int32(PLplot.Independent), Int32(PLplot.Default))

    for (i, xc) in enumerate(xlbls)
        if length(vlnc[i]) == 0
            scatter([xc], vlnw[i], '⚬')
        else
            xs = vcat(vlnw[i], -reverse(vlnw[i])) + xc
            ys = vcat(vlnc[i], reverse(vlnc[i]))
            PLplot.plline(length(xs), xs, ys)
        end
    end

end

function violin{T<:Real}(dct::Dict{T, Vector{T}}; trim=false, npoints=200, bw=0.0, width=0.8)
    xs = sort!(collect(keys(dct)))
    ys = [convert(Vector{PLplot.PLFLT}, dct[x]) for x in xs]
    violin(xs, ys, trim=trim, npoints=npoints, bw=bw, width=width)
end

function violin{T<:Real}(data::Matrix{T}; trim=false, npoints=200, bw=0.0, width=0.8)
    xs = collect(1:size(data,2))
    ys = [convert(Vector{PLplot.PLFLT}, data[:,i]) for i in 1:size(data,2)]
    violin(xs, ys, trim=trim, npoints=npoints, bw=bw, width=width)
end


# CI plot
#=
% ciplot(lower,upper)
% ciplot(lower,upper,x)
% ciplot(lower,upper,x,colour)
%
% Plots a shaded region on a graph between specified lower and upper confidence intervals (L and U).
% l and u must be vectors of the same length.
% Uses the 'fill' function, not 'area'. Therefore multiple shaded plots
% can be overlayed without a problem. Make them transparent for total visibility.
% x data can be specified, otherwise plots against index values.
=#
function ciplot{T<:Real}(x::AbstractVector{T}, data::AbstractMatrix{T}; wwidth=0.05, kwargs...)
    @assert size(data,1) == 3 "Data should have three rows: value, upper, lower"
    @assert size(data,2) == length(x) "Dimension mismatch"

    # wiskers
    ww = 0.05
    for i in 1:size(data,2)
        lw = data[2,i] # lower limit
        uw = data[3,i] # upper limit
        PLplot.pllsty(2)
        PLplot.pljoin(x[i], uw, x[i], lw)
        PLplot.pllsty(1)
        PLplot.pljoin(x[i]-wwidth, uw, x[i]+wwidth, uw)
        PLplot.pljoin(x[i]-wwidth, lw, x[i]+wwidth, lw)
    end
end
ciplot{T<:Real}(data::AbstractMatrix{T}; kwargs...) =
    ciplot(convert(Vector{Float64}, collect(1:size(data,2))), data; kwargs...)
ciplot!{T<:Real}(x::AbstractVector{T}, data::AbstractMatrix{T}; kwargs...) =
    ciplot(x, data; overlay=true, kwargs...)
ciplot!{T<:Real}(data::AbstractMatrix{T}; kwargs...) =
    ciplot(data; overlay=true, kwargs...)
