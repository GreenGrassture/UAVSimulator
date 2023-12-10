# For now these are all the same, but if I end up switching to 
# parameterizing by different A values this will make more sense

const A64 = UInt128(0xfc0072fa0b15f4fd)
const M = UInt128(2)^64

struct LCG
    A::UInt128
    privateSeed::UInt64
    M::UInt128
end

LCG() = LCG(A64, rand(UInt64), M)