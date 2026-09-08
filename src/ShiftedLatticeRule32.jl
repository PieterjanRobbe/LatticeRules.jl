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
```jldoctest; setup = :(using LatticeRules)
julia> lattice_rule = LatticeRule32([UInt32(1), UInt32(5)], 2, 8)
LatticeRule32{2}

julia> shifted_lattice_rule = ShiftedLatticeRule32(lattice_rule, [0.25, 0.5])
ShiftedLatticeRule32{2}

julia> getpoint(shifted_lattice_rule, 2)
2-element Vector{Float64}:
 0.5
 0.75

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
```jldoctest; setup = :(using LatticeRules)
julia> shifted_lattice_rule = ShiftedLatticeRule32(16)
ShiftedLatticeRule32{16}

julia> ndims(shifted_lattice_rule)
16

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
