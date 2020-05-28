# default generating vector files
const K_3600_32_file = joinpath(@__DIR__(), "..", "generating_vectors", "K_3600_32.txt")
const CKN_250_20_file = joinpath(@__DIR__(), "..", "generating_vectors", "CKN_250_20.txt")

# re-compile if these files change
include_dependency(K_3600_32_file)
include_dependency(CKN_250_20_file)

# default generating vectors
const K_3600_32 = read32(K_3600_32_file)
const CKN_250_20 = read32(CKN_250_20_file)

