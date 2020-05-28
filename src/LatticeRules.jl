module LatticeRules

import Random

export AbstractLatticeRule, LatticeRule, LatticeRule32, ShiftedLatticeRule, ShiftedLatticeRule32, getpoint, unsafe_getpoint, unsafe_getpoint!, ndims, length

export K_3600_32, K_3600_32_file, CKN_250_20, CKN_250_20_file # from lattice_data.jl

for file in ["common", "lattice_data", "LatticeRule32", "ShiftedLatticeRule32"]
    include(string(file, ".jl"))
end

end # module
