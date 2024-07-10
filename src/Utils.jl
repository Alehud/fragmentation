export lshift, rshift, arr2int, brint, flip_bit, get_bit, vertex_labels, remove_vacuum!, arrow_func, int2label, int2str, mymod, plaquette_idx, 
row_idx, column_idx, draw_grid, issubarray


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

Convert a boolean vector `a` of length `n` to an integer.

# Arguments
- `a::Vector{Bool}`: boolean vector

# Returns
- `Integer`: an integer with the binary representation stored in the vector `a`
"""
function arr2int(a::Vector{<:Integer})::Integer
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
function vertex_labels(states::Vector{<:Vector{<:Integer}}, legend::Vector{<:Pair{<:Integer, String}} = nothing)
    if isnothing(legend)
        return [join(string.(state)) for state in states]
    else
        return [join(replace(state, legend...)) for state in states]
    end
end


function remove_vacuum!(state_labels::Vector{String}, g::Graph)
    i = 1
    while i ≤ length(state_labels)
        vs = findall(x -> replace(x, "_" => "") == replace(state_labels[i], "_" => ""), state_labels)
        merge_vertices!(g, vs)
        deleteat!(state_labels, vs[2:end])
        state_labels[i] = replace(state_labels[i], "_" => "")
        i += 1
    end
end


function arrow_func(s, d, y)
    if abs(y[s]-y[d]) == 1 || abs(s-d) == 1
        return [0, 0]
    else
        return [20, 25]
    end
end


# function serialize_py(foldername::AbstractString, filename::AbstractString, value)
#     println("Serializing $(filename)...")
#     flush(stdout)
#     serialize("$(foldername)/$(filename)", value)
#     println("   $(filename) serialized")
#     println("Converting $(filename) to python...")
#     flush(stdout)
#     if typeof(value) <: Dict
#         @pywith pybuiltin("open")("$(foldername)/python_data/$(filename[1:end-4]).pkl","wb") as f begin
#             pickle.dump(PyDict(Dict(Tuple.(keys(value)) .=> values(value))), f)
#         end
#     elseif typeof(value) <: Vector
#         @pywith pybuiltin("open")("$(foldername)/python_data/$(filename[1:end-4]).pkl","wb") as f begin
#             pickle.dump(PyObject(value), f)
#         end
#     end
#     println("   $(filename) pickled")
#     flush(stdout)
# end


function int2label(t::Integer)
    if t == 0
        t_str = "0"
    else
        ex = Int(floor(log10(t)))
        if ex ≤ 2
            t_str = string(t)
        else
            base = round(t / 10^(ex-2)) / 10^2
            if base == 1.0
                t_str = "10^{$(ex)}"
            else
                if isinteger(base)
                    base = Int(base)
                end
                t_str = "$(base) \\cdot 10^{$(ex)}"
            end
        end
    end
    return t_str
end

function int2str(t::Integer)
    if t == 0
        t_str = "0"
    else
        ex = Int(floor(log10(t)))
        if ex ≤ 2
            t_str = string(t)
        else
            base = round(t / 10^(ex-2)) / 10^2
            if base == 1.0
                t_str = "10^$(ex)"
            else
                if isinteger(base)
                    base = Int(base)
                end
                t_str = "$(base)p10^$(ex)"
            end
        end
    end
    return t_str
end


function mymod(a::Integer, n::Integer)
    return mod(a-1, n) + 1
end

function mymod(v::Vector{<:Integer}, n::Integer)
    return (a -> mymod(a, n)).(v)
end

function plaquette_idx(column::Integer, row::Integer; Lx::Integer, Ly::Integer)
    return Int64[2((row-1)*Lx + column-1) + 1, 2((row-1)*Lx + column-1) + 2, 2((mymod(row+1, Ly)-1)*Lx + column-1) + 1, 2((row-1)*Lx + mymod(column+1,Lx)-1) + 2]
end

# function star_idx(column::Integer, row::Integer; Lx::Integer, Ly::Integer)
#     return Int64[2((row-1)*Lx + column-1) + 1, 2((row-1)*Lx + column-1) + 2, 2((mymod(row+1, Ly)-1)*Lx + column-1) + 1, 2((row-1)*Lx + mymod(column+1,Lx)-1) + 2]
# end

function row_idx(row::Integer; Lx::Integer)
    return Int64[2(row-1)*Lx + 1:2:2(row*Lx -1) + 1;]
end

function column_idx(column::Integer; Lx::Integer, Ly::Integer)
    return Int64[2column:2Lx:2((Ly-1)*Lx + column-1) + 2;]
end

function draw_grid(state::Vector{Int8}; Lx::Integer, Ly::Integer, legend::Union{Dict, Nothing}=nothing)
    if !isnothing(legend)
        s = replace(state, legend...)
    else
        s = state
    end
    for y in 0:Ly-1
        for x in 0:Lx-1
            print("+--$(s[2x+1 + 2y*Lx])--")
        end
        println("+")
        for x in 1:Lx
            print("|     ")
        end
        println("|")
        for x in 0:Lx-1
            print("$(s[2x+2 + 2y*Lx])     ")
        end
        println("$(s[2 + 2y*Lx])")
        for x in 1:Lx
            print("|     ")
        end
        println("|")
    end
    for x in 0:Lx-1
        print("+--$(s[2x+1])--")
    end
    println("+")
end


function issubarray(needle, haystack)
    getView(vec::Vector{Int8}, i::Integer, len::Integer) = view(vec, i:i+len-1)
    ithview(i::Integer) = getView(haystack, i, length(needle))
    return any(i -> ithview(i) == needle, 1:length(haystack)-length(needle)+1)
end