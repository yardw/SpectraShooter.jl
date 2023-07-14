module tmpmod
export Spectra

"""
    Spectra

This is a simple wrapper around a 2D array containing the approximate position of multiple spectra lines. The first dimension is by default the eigenvalue(like the masses of particles), while the second is by default the parameter on searching(e.g. one of the free parameters in the model), so that it shows how the parameter would affect the spectra.
Components:
- `data`: the 2D array containing the approximate position of multiple spectra lines.
- `lograngeeigval`: the range of the eigenvalue, in log scale.
- `lograngeparam`: the range of the parameter, in log scale.
- `resolusioneigval`: the resolution of the eigenvalue, i.e. the number of points in the first dimension.
- `resolusionparam`: the resolution of the parameter, i.e. the number of points in the second dimension.
"""
    struct Spectra <: AbstractArray{Float64,2}
        data::Array{Float64,2}
        # keep track of the range of the eigenvalue and the parameter
        lograngeeigval::Tuple{Float64,Float64}
        lograngeparam::Tuple{Float64,Float64}
        # constructor
        function Spectra(rangeeigval::Tuple{Float64,Float64}, 
                         rangeparam::Tuple{Float64,Float64}, 
                         (resolusioneigval=64, resolusionparam=64)::Tuple{Int64,Int64}, 
                         logscaled::Bool=true )
            new(zeros((resolusioneigval, resolusionparam)), 
                rangeeigval, rangeparam)
        end
    end
    # overload of the `:AbstractArray` essential interfaces
    Base.size(s::Spectra) = size(s.data)
    Base.IndexStyle(::Type{<:Spectra}) = IndexCartesian()
    Base.getindex(s::Spectra, i::Int, j::Int) = s.data[i,j]
    Base.setindex!(s::Spectra, v, i::Int, j::Int) = s.data[i,j] = v
end