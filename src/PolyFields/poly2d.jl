abstract type AbstractBivariatePolynomial end

struct BivariatePolynomial{F,N,X,M,Y} <: AbstractBivariatePolynomial
    p::ImmutablePolynomial{ImmutablePolynomial{F,X,N},Y,M}
    function BivariatePolynomial{F}(t::Tuple) where F
            all(typeof(u)<:Tuple for u in t) || throw(ArgumentError("A tuple of tuples is expected"))
            N = maximum(length.(t))
            M = length(t)
            q = Tuple(ImmutablePolynomial(_zerofill(F.(u),N),:x) for u in t)
            p = ImmutablePolynomial(q,:y)
            new{F,N,:x,M,:y}(p)
    end
end

_zerofill(t::NTuple{K,T},N) where {K,T} = K>=N ? t : (t...,zeros(T,N-K)...)

struct TensorPolynomial{F,X,N,Y,M} <: AbstractBivariatePolynomial
    px::ImmutablePolynomial{F,X,N}
    py::ImmutablePolynomial{F,Y,M}
end

