abstract type AbstractBivariatePolynomial <: AbstractPolynomial end

struct BivariatePolynomial{F,N,X,M,Y} <: AbstractPolynomial
    p::ImmutablePolynomial{ImmutablePolynomial{F,X,N},Y,M}
end

_zerofill(t::NTuple{K,T},N) where {K,T} = K>=N ? t : (t...,zeros(T,N-K)...) 