using Clustering
using LMCLUS
using NetCDF
using JSON

include("transformations.jl")

function labelManifolds(Ms::Array{Manifold}, s::Int = sum(map(outdim, Ms)))
    C = zeros(Int, s)
    for i = 1:length(Ms)
        idxs = labels(Ms[i])
        for j = 1:length(idxs)
            idx = idxs[j]
            try
                if C[idx] != 0
                    warn("Assigning new index $(i) to point $(idxs[j]) that has already index $(C[idxs[j]])")
                end
                C[idxs[j]] = i
            catch
                error("Index $(idx) out of bounds 1:$(s)")
            end
        end
    end
    return C
end

function labelManifolds(Ms::Array{Manifold}, extidxs::Vector{Int}, s::Int = sum(map(outdim, Ms)))
    C = zeros(Int, s)
    for i = 1:length(Ms)
        idxs = labels(Ms[i])
        for j = 1:length(idxs)
            idx = extidxs[idxs[j]]
            try
                if C[idx] != 0
                    warn("Assigning new index $(i) to point $(idxs[j]) that has already index $(C[idxs[j]])")
                end
                C[extidxs[idxs[j]]] = i
            catch
                error("Index $(idx) out of bounds 1:$(s)")
            end
        end
    end
    return C
end

function loadClusters(::Type{Manifold}, fn)
    fid = NetCDF.open(fn) #  readdimvar=true
    μs = NetCDF.readvar(fid, "mu")
    labels = NetCDF.readvar(fid, "labels")
    N, m = size(μs)
    L = 1:length(labels)
    manifolds = Manifold[]
    for i in 1:m
        mname = "M$i"
        μ = μs[:,i]
        lbls = L[labels .== i]
        separation = Separation()
        try
            projection = NetCDF.readvar(fid, mname)
            dim = size(projection,2)
            separation.depth = fid.vars[mname].atts["depth"]
            separation.discriminability = fid.vars[mname].atts["discriminability"]
            separation.threshold = fid.vars[mname].atts["θ"]
            push!(manifolds, Manifold(dim,μ,projection,lbls,separation))
        catch
            dim = 0
            projection = ones(dim,dim)
            push!(manifolds, Manifold(dim,μ,projection,lbls,separation))
        end
    end
    NetCDF.close(fid)
    return manifolds
end

function loadClusters(::Type{KmeansResult}, fn)
    centers = ncread(fn, "mu")
    labels = ncread(fn, "labels")
    return centers, labels
end

function saveClusters(fn, clusters::KmeansResult, t::DataTransform.UnitRangeTransform; params...)
    params = Dict{Symbol,Any}(params)
    cdims = size(clusters.centers)
    filename = "$(fn).nc"
    indexes = clusters.assignments
    d1 = NcDim("N",cdims[1];atts=Dict{Any,Any}("longname"=>"Dimension of the dataset"))
    d2 = NcDim("k",cdims[2];atts=Dict{Any,Any}("longname"=>"Number of clusters"))
    d3 = NcDim("n",length(indexes);atts=Dict{Any,Any}("longname"=>"Size of the dataset"))

    tscale = NcVar("t_scale",d1)
    tmin = NcVar("t_min",d1)
    μ = NcVar("μ",[d1,d2])
    lbls = NcVar("labels",d3, t=Int32)

    nc = NetCDF.create(filename,[tscale,tmin,μ,lbls],gatts=Dict{Any,Any}(
            "rng_seed" => params[:rng_seed],
            "k"  => cdims[2],
            "lon" => params[:res][1],
            "lat" => params[:res][2],
            "t_unit" => int(t.unit)
        )
    )
    NetCDF.putvar(nc,"t_scale", t.scale)
    NetCDF.putvar(nc,"t_min", t.min)
    NetCDF.putvar(nc,"μ",clusters.centers)
    NetCDF.putvar(nc,"labels",indexes)
    NetCDF.close(nc)
end

function label_manifolds(Ms::Array{Manifold}, X::AbstractArray)
    C = zeros(Int, size(X,2))
    for i = 1:length(Ms)
        C[labels(Ms[i])] = i
    end
    return C
end

function manifold_var(mname::AbstractString, m::Manifold, d1::NcDim)
    d2 = NcDim(mname,indim(m);atts=Dict{Any,Any}("longname"=>"Dimension of the manifold"))
    return indim(m) == 0 ? NcVar(mname,[d1,d2], atts=Dict{Any,Any}("longname"=>"Noise")) :
        NcVar(mname,[d1,d2]; atts=Dict{Any,Any}(
        "θ" => m.separation.threshold,
        "depth" => m.separation.depth,
        "discriminability" => m.separation.discriminability )
    )
end

function saveClusters(fn, manifolds::Array{Manifold}, p::LMCLUSParameters,
                        L::Vector{Int}, t::DataTransform.UnitRangeTransform; gparams...)
    filename = "$(fn).nc"

    # Form global attributes
    gparams = Dict{Symbol,Any}(gparams)
    params = Dict{Any,Any}(JSON.Parser.parse(JSON.json(p)))
    params["manifolds"] = length(manifolds)
    for (k,v) in params
        if isa(v,Bool)
            params[k] = int(v)
        end
    end
    params["lon"] = gparams[:res][1]
    params["lat"] = gparams[:res][2]
    params["t_unit"] = int(t.unit)

    #n = sum(map(x->outdim(x), manifolds))
    N = size(projection(manifolds[1]),1)

    d1 = NcDim("N",N;atts=Dict{Any,Any}("longname"=>"Dimension of the dataset"))
    d2 = NcDim("n",length(L);atts=Dict{Any,Any}("longname"=>"Size of the dataset"))
    d3 = NcDim("m",length(manifolds);atts=Dict{Any,Any}("longname"=>"Number of LM clusters"))

    tscale = NcVar("t_scale",d1)
    tmin = NcVar("t_min",d1)
    μ = NcVar("μs",[d1,d3])
    lbls = NcVar("labels",d2, t=Int32)
    mvars = Array(NcVar, length(manifolds))
    for i in 1:length(manifolds)
        mvars[i] = manifold_var("M$i", manifolds[i], d1)
    end

    nc = NetCDF.create(filename,[tscale,tmin,μ,lbls,mvars],gatts=params)
    NetCDF.putvar(nc,"t_scale", t.scale)
    NetCDF.putvar(nc,"t_min", t.min)
    NetCDF.putvar(nc,"μs",hcat([mean(m) for m in manifolds]...))
    NetCDF.putvar(nc,"labels", L)
    for i in 1:length(manifolds)
        if indim(manifolds[i]) > 0
            NetCDF.putvar(nc,"M$i",projection(manifolds[i]))
        end
    end
    NetCDF.close(nc)
end

function loadTransform(fn)
    s = ncread(fn, "t_scale")
    m = ncread(fn, "t_min")
    return DataTransform.UnitRangeTransform(length(s), true, m, s)
end

function loadGlobal(fn)
    fid = NetCDF.open(fn, readdimvar=true)
    attr = Dict{String,Any}(fid.gatts)
    NetCDF.close(fid)
    return attr
end

function contingency(L1::Vector, L2::Vector)
    @assert length(L1) == length(L2)
    c1 = length(unique(L1)) # classes
    c2 = length(unique(L2)) # clusters
    ct = zeros(Int, c1, c2)
    for i in 1:length(L1)
        ct[L1[i], L2[i]] += 1
    end
    return ct
end

# Accepts contingency table with row as classes, c, and columns as clusters, k.
function V_measure(A::Matrix; β = 1.0)
    C, K = size(A)
    N = sum(A)

    # Homogeneity
    hck = 0.0
    for k in 1:K
        d = sum(A[:,k])
        for c in 1:C
            if A[c,k] != 0 && d != 0
                hck += log(A[c,k]/d) * A[c,k]/N
            end
        end
    end
    hck = -hck

    hc = 0.0
    for c in 1:C
        n = sum(A[c,:]) / N
        if n != 0.0
            hc += log(n) * n
        end
    end
    hc = -hc

    h = (hc == 0.0 || hck == 0.0) ? 1 : 1 - hck/hc

    # Completeness
    hkc = 0.0
    for c in 1:C
        d = sum(A[c,:])
        for k in 1:K
            if A[c,k] != 0 && d != 0
                hkc += log(A[c,k]/d) * A[c,k]/N
            end
        end
    end
    hkc = -hkc

    hk = 0.0
    for k in 1:K
        n = sum(A[:,k]) / N
        if n != 0.0
            hc += log(n) * n
        end
    end
    hk = -hk

    c = (hk == 0.0 || hkc == 0.0) ? 1 : 1 - hkc/hk

    # V-measure
    V_β = (1 + β)*h*c/(β*h + c)
    return V_β
end

function count_clust(lbls::Vector{Int32})
    l = length(lbls)
    c = zeros(unique(lbls))
    zidx = 0 in unique(lbls) ? 1 : 0

    for i in 1:l
        c[lbls[i]+zidx] += 1
    end
    z = zidx == 1 ? shift!(c) : 0
    return c, z
end