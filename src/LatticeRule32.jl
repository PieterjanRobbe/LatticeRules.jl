#
# LatticeRule32
#
struct LatticeRule32{s} <: AbstractLatticeRule{s}
    z::Vector{UInt32} # generating vector
    n::Int64 # max number of points in the lattice rule
end

# default lattice rule type
const LatticeRule = LatticeRule32

# access max number of points in the lattice
Base.length(lattice_rule::LatticeRule32) = lattice_rule.n

# uinttype
uinttype(::LatticeRule32) = UInt32

"""
    LatticeRule32(z, s, n)
    LatticeRule32(z, s)
    LatticeRule32(z)

Returns a rank-1 lattice rule in `s` dimensions with generating vector `z` and at most `n` points.

When no maximum number of points `n` is provided, we assume `n = typemax(UInt32) = 2^32 - 1`. When no number of dimensions `s` is provided, we assume `s = length(z)`. 

!!! info

    Technically, we return an extensible lattice sequence where the `k`-th point is transformed using the gray coded radical inverse function. This has the advantage that we can add points to the lattice without changing the already computed points.

More generating vectors can be found online [here](https://web.maths.unsw.edu.au/~fkuo/lattice/index.html) or [here](https://people.cs.kuleuven.be/~dirk.nuyens/qmc-generators/).

# Examples
```jldoctest; setup = :(using LatticeRules; import Random; Random.seed!(1))
julia> lattice_rule = LatticeRule32([UInt32(1), UInt32(5)], 2, 8) # Fibonacci lattice
LatticeRule32{2}

julia> getpoint(lattice_rule, 2)
2-element Array{Float64,1}:
 0.25
 0.25

```
See also: [`getpoint`](@ref), [`ShiftedLatticeRule32`](@ref)
"""
LatticeRule32(z::Vector{UInt32}) = LatticeRule32(z, length(z)) # specify generating vector

# specify generating vector and number of dimensions
LatticeRule32(z::Vector{UInt32}, s::Integer) = LatticeRule32(z, s, typemax(UInt32) + 1)

# specify generating vector, number of dimensions and maximum number of points
function LatticeRule32(z::Vector{UInt32}, s::Integer, n::Integer)
    s > 0 || throw(ArgumentError("number of dimensions s must be larger than 0"))
    s ≤ length(z) || throw(ArgumentError("number of dimensions s must be less than or equal to the length of the generating vector z"))
    n > 0 || throw(ArgumentError("maximum number of points n must be larger than 0"))
    n ≤ typemax(UInt32) + 1 || throw(ArgumentError("maximum number of points n must be less than or equal to 2^32, consider implementing a LatticeRule64 type"))
    LatticeRule32{s}(view(z, 1:s), n)
end

"""
    LatticeRule32(file, s, n)
    LatticeRule32(file, s)
    LatticeRule32(file)

Returns a rank-1 lattice rule in `s` dimensions with generating vector `z` read from the file `file` and with at most `n` points.

# Examples
```jldoctest; setup = :(using LatticeRules)
julia> z_file = K_3600_32_file;

julia> lattice_rule = LatticeRule32(z_file, 16)
LatticeRule32{16}

julia> getpoint(lattice_rule, 123)
16-element Array{Float64,1}:
 0.8671875
 0.9609375
 0.6015625
 0.8984375
 0.6484375
 0.6328125
 0.3203125
 0.2890625
 0.0234375
 0.1015625
 0.7890625
 0.0703125
 0.6953125
 0.0234375
 0.1171875
 0.0859375

```
See also: [`getpoint`](@ref), [`ShiftedLatticeRule32`](@ref)
"""
LatticeRule32(file::AbstractString) = LatticeRule32(read32(file)) # specify file containing generating vector

# specify file containting generating vector and number of dimensions
LatticeRule32(file::AbstractString, s::Integer) = LatticeRule32(read32(file), s)

# specify file containing generating vector, number of dimensions and maximum number of points
LatticeRule32(file::AbstractString, s::Integer, n::Integer) = LatticeRule32(read32(file), s, n)

"""
    LatticeRule32(s)

Returns a rank-1 lattice rule in `s` dimensions that uses a default generating vector with order-2 weights.

# Examples
```jldoctest; setup = :(using LatticeRules)
julia> lattice_rule = LatticeRule32(16)
LatticeRule32{16}

julia> getpoint(lattice_rule, 123)
16-element Array{Float64,1}:
 0.8671875
 0.5390625
 0.6015625
 0.3671875
 0.6796875
 0.8203125
 0.3046875
 0.8515625
 0.7109375
 0.6328125
 0.5703125
 0.2578125
 0.6953125
 0.0390625
 0.2421875
 0.4453125

```
See also: [`getpoint`](@ref), [`ShiftedLatticeRule32`](@ref)
"""
function LatticeRule32(s::Integer) # specify number of dimensions only
    s ≤ 3600 || throw(ArgumentError("number of dimensions s must be less than or equal to 3600, please supply your own generating vector z"))
    s ≤ 250 ? LatticeRule32(CKN_250_20, s, 2^20) : LatticeRule32(K_3600_32, s)
end

# in-place version of unsafe_getpoint (with 0 memory allocations)
@inline function unsafe_getpoint!(x::Vector{<:AbstractFloat}, lattice_rule::LatticeRule32, k::UInt32)
    ϕ_k = reversebits(k) * 2.0^(-32) # gray coded radical inverse function in base 2
    @inbounds for i in 1:length(x)
        x[i] = ϕ_k * lattice_rule.z[i]
        x[i] -= floor(x[i]) # mod 1
    end
    x
end

# fancy printing
Base.show(io::IO, lattice_rule::LatticeRule32{s}) where s = print(io, string("LatticeRule32{", s, "}"))
