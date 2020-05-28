#
# AbstractLatticeRule
#
abstract type AbstractLatticeRule{s} <: Random.AbstractRNG end

# number of dimensions
Base.ndims(::AbstractLatticeRule{s}) where s = s::Int

# size
Base.size(lattice_rule::AbstractLatticeRule) = (length(lattice_rule), )

# reverse bits (https://graphics.stanford.edu/~seander/bithacks.html#ReverseParallel)
reversebits(u::UInt32) = begin
    u = ((u >> 1) & 0x55555555) | ((u & 0x55555555) << 1)
    u = ((u >> 2) & 0x33333333) | ((u & 0x33333333) << 2)
    u = ((u >> 4) & 0x0F0F0F0F) | ((u & 0x0F0F0F0F) << 4)
    u = ((u >> 8) & 0x00FF00FF) | ((u & 0x00FF00FF) << 8)
    u = ( u >> 16             ) | ( u               << 16)
end

# read contents of file containing a generating vector and convert it to a Vector of UInt32's
read32(file::AbstractString) = parse.(UInt32, readlines(file))

"""
    getpoint(lattice_rule, k)

Get the `k`-th point of the lattice rule `lattice_rule`.

!!! note

    An alternative syntax is `getindex(lattice_rule, k)` or `lattice_rule[k]`, this allows you to write the one-liner `Q = mean(f.(lattice_rule[0:N-1]))` for the quasi-Monte Carlo estimator for ``E[f]``.

```jldoctest; setup = :(using LatticeRules)
julia> lattice_rule = LatticeRule32(2)
LatticeRule32{2}

julia> getpoint(lattice_rule, 3)
2-element Array{Float64,1}:
 0.75
 0.25

```
See also: [`LatticeRule32`](@ref), [`ShiftedLatticeRule32`](@ref)
"""
@inline function getpoint(lattice_rule::AbstractLatticeRule, k::Number) # get the k-th point of the lattice sequence
    0 ≤ k < length(lattice_rule) || throw(BoundsError(lattice_rule, k))
    unsafe_getpoint(lattice_rule, convert(uinttype(lattice_rule), k))
end

# get the k-th point of the sequence without bounds checking for 32 bit integers
@inline unsafe_getpoint(lattice_rule::AbstractLatticeRule{s}, k::UInt32) where s = begin
    x = Vector{Float64}(undef, s)
    unsafe_getpoint!(x, lattice_rule, k) # dispatch to AbstractLatticeRule subtype
end

# make LatticeRule iterable
Base.iterate(lattice_rule::AbstractLatticeRule, state=uinttype(lattice_rule)(0)) = state ≥ length(lattice_rule) ? nothing : (getpoint(lattice_rule, state), state + uinttype(lattice_rule)(1))
Base.eltype(::Type{<:AbstractLatticeRule}) = Vector{Float64}

# enable lattice_rule[i] access
Base.getindex(lattice_rule::AbstractLatticeRule, i::Number) = getpoint(lattice_rule, i)
Base.getindex(lattice_rule::AbstractLatticeRule, I) = [lattice_rule[i] for i in I]
Base.firstindex(lattice_rule::AbstractLatticeRule) = uinttype(lattice_rule)(0)
Base.lastindex(lattice_rule::AbstractLatticeRule) = uinttype(lattice_rule)(length(lattice_rule) - 1)
