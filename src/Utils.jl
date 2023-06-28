export lshift, rshift, arr2int, brint, flip_bit, get_bit, vertex_labels


"""
    lshift(x, d, n)

Circular left bit shift of number `x` by `d` bits with the total length of `n`.

# Arguments
- `x::Integer`: an integer number
- `d::Integer`: by how many bits to shift
- `n::Integer`: the total length of the number in bits

# Returns
- `Integer`: circularly shifted by d bits to the left
"""
function lshift(x::Integer, d::Integer, n::Integer)::Integer
    return ((x << d) % (1 << n)) | (x >> (n - d))
end

function lshift(x::BigInt, d::Integer, n::Integer)::Integer
    return ((x << d) % (BigInt(1) << n)) | (x >> (n - d))
end



"""
    rshift(x, d, n)

Circular right bit shift of the integer `x` by `d` bits with the total length of `n`.

# Arguments
- `x::Integer`: an integer number
- `d::Integer`: by how many bits to shift
- `n::Integer`: the total length of the number in bits

# Returns
- `Integer`: circularly shifted by `d` bits to the right
"""
function rshift(x::Integer, d::Integer, n::Integer)::Integer
    return (x >> d) | (x << (n - d)) % (1 << n)
end

function rshift(x::BigInt, d::Integer, n::Integer)::Integer
    return (x >> d) | (x << (n - d)) % (BigInt(1) << n)
end


"""
    arr_to_int(a)

Convert numpy boolean vector `a` of length `n` to an integer.

# Arguments
- `a::Vector{Bool}`: boolean vector

# Returns
- `Integer`: an integer with the binary representation stored in the vector `a`
"""
function arr2int(a::Vector{Bool})::Integer
    if length(a) > 63
        return sum([p*2^BigInt(i-1) for (i, p) in enumerate(a)])
    else
        return sum([p*2^(i-1) for (i, p) in enumerate(a)])
    end
end


"""
    brint(x, n)

Print integers in binary (with length `n`)

# Arguments
- `x::Union{Integer, Vector{<:Integer}}`: an integer number or a Vector of integers
- `n::Integer=8`: the length in bits
"""
function brint(x::Union{Integer, Vector{<:Integer}}, n::Integer=8)::Nothing
    if isa(x, Integer)
        println(digits(Bool, x, base=2, pad=n))
    elseif isa(x, Vector{<:Integer})
        for y in x
            println(digits(Bool, y, base=2, pad=n))
        end
    end
end


"""
    flip_bit(x, k)

Flip `k`'th bit of integer `x`.

# Arguments
- `x::Integer`: an integer number
- `k::Integer`: bit to flip (starts from 0)

# Returns
- `Integer`: integer `x` with `k`'th bit flipped
"""
function flip_bit(x::Integer, k::Integer)::Integer
    return x ⊻ (1 << k)
end

function flip_bit(x::BigInt, k::Integer)::Integer
    return x ⊻ (BigInt(1) << k)
end


"""
    get_bit(x, k)

Get k'th bit of integer x.

# Arguments
- `x::Integer`: an integer number
- `k::Integer`: bit to get (starts from 0)

# Returns
- `Bool`: `k`'th bit of integer `x`
"""
function get_bit(x::Integer, k::Integer)::Bool
    return (x >> k) & 1
end

function get_bit(x::BigInt, k::Integer)::Bool
    return (x >> k) & BigInt(1)
end


"""
    vertex_labels(states, legend)

Create string labels out of the vector with states.

# Arguments
- `states::Vector{Vector{<:Integer}}`: vector with states
- `legend::Vector{Pair{String, String}} = nothing`: rules for substituting integers with symbols

# Returns
- `Vector{String}`: string labels
"""
function vertex_labels(states::Vector{<:Vector{<:Integer}}, legend::Vector{Pair{String, String}} = nothing)
    if isnothing(legend)
        return [join(string.(state)) for state in states]
    else
        return [replace(join(string.(state)), legend...) for state in states]
    end
end
