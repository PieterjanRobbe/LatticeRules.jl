#
# ShiftedLatticeRule32
#
struct ShiftedLatticeRule32{s, L, V} <: AbstractLatticeRule{s}
    lattice_rule::L
    Δ::V
end

# default shifted lattice rule type
const ShiftedLatticeRule = ShiftedLatticeRule32

# access max number of points in the lattice
Base.length(shifted_lattice_rule::ShiftedLatticeRule32) = length(shifted_lattice_rule.lattice_rule)

# uinttype
uinttype(::ShiftedLatticeRule32) = UInt32

"""
    ShiftedLatticeRule32(lattice_rule)
    ShiftedLatticeRule32(lattice_rule, shift)

Returns a shifted rank-1 lattice rule based on the lattice rule `lattice_rule` using the random shift `shift`. If no random shift is provided, we use `shift = rand(length(lattice_rule))`.

# Examples
```jldoctest; setup = :(using LatticeRules; import Random; Random.seed!(1))
julia> lattice_rule = LatticeRule32(16)
LatticeRule32{16}

julia> shifted_lattice_rule = ShiftedLatticeRule32(lattice_rule)
ShiftedLatticeRule32{16}

julia> getpoint(shifted_lattice_rule, 0)
16-element Array{Float64,1}:
 0.23603334566204692
 0.34651701419196046
 0.3127069683360675
 0.00790928339056074
 0.4886128300795012
 0.21096820215853596
 0.951916339835734
 0.9999046588986136
 0.25166218303197185
 0.9866663668987996
 0.5557510873245723
 0.43710797460962514
 0.42471785049513144
 0.773223048457377
 0.2811902322857298
 0.20947237319807077

```
See also: [`LatticeRule32`](@ref), [`getpoint`](@ref)
"""
ShiftedLatticeRule32(lattice_rule::LatticeRule32{s}) where s = ShiftedLatticeRule32(lattice_rule, rand(s)) # specify lattice rule

# specify lattice rule and random shift
function ShiftedLatticeRule32(lattice_rule::LatticeRule32{s}, Δ::Vector{<:AbstractFloat}) where s
    length(Δ) == ndims(lattice_rule) || throw(DimensionMismatch("length of the random shift vector must be equal to the number of dimensions of the lattice rule, expected $(ndims(lattice_rule)), got $(length(Δ))"))
    all(0 .≤ Δ .≤ 1) || throw(ArgumentError("random shift vector must contain uniformly distributed random numbers"))
    ShiftedLatticeRule32{s, typeof(lattice_rule), typeof(Δ)}(lattice_rule, Δ)
end

"""
    ShiftedLatticeRule32(s)

Returns a shifted rank-1 lattice rule in `s` dimensions that uses a default generating vector with order-2 weights and a randomly generated shift vector.

# Examples
```jldoctest; setup = :(using LatticeRules; import Random; Random.seed!(1))
julia> shifted_lattice_rule = ShiftedLatticeRule32(16)
ShiftedLatticeRule32{16}

julia> shifted_lattice_rule[0]
16-element Array{Float64,1}:
 0.23603334566204692
 0.34651701419196046
 0.3127069683360675
 0.00790928339056074
 0.4886128300795012
 0.21096820215853596
 0.951916339835734
 0.9999046588986136
 0.25166218303197185
 0.9866663668987996
 0.5557510873245723
 0.43710797460962514
 0.42471785049513144
 0.773223048457377
 0.2811902322857298
 0.20947237319807077

```
See also: [`getpoint`](@ref), [`ShiftedLatticeRule32`](@ref)
"""
ShiftedLatticeRule32(s::Integer) = ShiftedLatticeRule32(LatticeRule32(s)) # specify number of dimensions only

# in-place version of unsafe_getpoint (with 0 memory allocations)
@inline function unsafe_getpoint!(x::Vector{<:AbstractFloat}, shifted_lattice_rule::ShiftedLatticeRule32, k::UInt32)
    ϕ_k = reversebits(k) * 2.0^(-32) # gray coded radical inverse function in base 2
    @inbounds for i in 1:length(x)
        x[i] = ϕ_k * shifted_lattice_rule.lattice_rule.z[i] + shifted_lattice_rule.Δ[i]
        x[i] -= floor(x[i]) # mod 1
    end
    x
end

# fancy printing
Base.show(io::IO, shifted_lattice_rule::ShiftedLatticeRule32{s}) where s = print(io, string("ShiftedLatticeRule32{", s, "}"))
