
struct ConicObjVAF{T}
    mapping::OrderedDict{UInt64, MOI.VectorAffineFunction{T}}
end

function make_VAF(unique_conic_forms::UniqueConicForms{T}, id::UInt64, val::Tuple{Value, Value}) where T
    re, im = val
    sz = size(re, 1)
    @assert sz == size(im, 1)
    terms = MOI.VectorAffineTerm{T}[]
    b = spzeros(T, sz)
    var = unique_conic_forms.id_to_variables[id]
    terms = MOI.VectorAffineTerm{T}[]
    b = spzeros(T, sz)
    if id == objectid(:constant)
        for l in 1:sz
            b[l] = ifelse(val[1][l] == 0, val[2][l],  val[1][l])
        end
    else
            var_indices = unique_conic_forms.id_to_indices[id]
            if unique_conic_forms.id_to_variables[id].sign == ComplexSign()
                l = length(var_indices) ÷ 2
                # We create MOI terms and add them to `terms`.
                # Real part:
                add_terms!(terms, val[1], 1 :  sz, var_indices[1:l])

                # Imaginary part:
                add_terms!(terms, val[2], 1 :  sz, var_indices[l+1 : end])

            else
                # Only a real part:
                add_terms!(terms, val[1], 1 :  sz, var_indices)
            end
        end

    return MOI.VectorAffineFunction{T}(terms, b)
end

function make_ConicObjVAF!(unique_conic_forms::UniqueConicForms{T}, o::ConicObj) where T
    d = OrderedDict( id => make_VAF(unique_conic_forms, id, val) for (id, val) in pairs(o)  )
    return  ConicObjVAF{T}(d)
end

Base.iterate(c::ConicObjVAF, s...) = iterate(c.mapping, s...)
Base.keys(c::ConicObjVAF) = keys(c.mapping)
Base.haskey(c::ConicObjVAF, var::Integer) = haskey(c.mapping, UInt64(var))
Base.getindex(c::ConicObjVAF, var::Integer) = c.mapping[UInt64(var)]
Base.setindex!(c::ConicObjVAF, val, var::Integer) = setindex!(c.mapping, val, UInt64(var))
Base.copy(c::ConicObjVAF) = ConicObjVAF(copy(c.mapping))



# # helper function to negate conic objectives
# # works by changing each (key, val) pair to (key, -val)
# function Base.:-(c::ConicObjVAF)
#     new_obj = copy(c)
#     for var in keys(new_obj)
#         x1 = new_obj[var][1]*(-1)
#         x2 = new_obj[var][2]*(-1)
#         new_obj[var] = (x1,x2)
#     end
#     return new_obj
# end

# function Base.:+(c::ConicObj, d::ConicObj)
#     new_obj = copy(c)
#     for var in keys(d)
#         if !haskey(new_obj, var)
#             new_obj[var] = d[var]
#         else
#             # .+ does not preserve sparsity
#             # need to override behavior
#             if size(new_obj[var][1]) == size(d[var][1])
#                 x1 = new_obj[var][1] + d[var][1]
#             else
#                 x1 = broadcast(+, new_obj[var][1], d[var][1])
#             end
#             if size(new_obj[var][2]) == size(d[var][2])
#                 x2 = new_obj[var][2] + d[var][2]
#             else
#                 x2 = broadcast(+, new_obj[var][2], d[var][2])
#             end
#             new_obj[var] = (x1,x2)
#         end
#     end
#     return new_obj
# end

# function get_row(c::ConicObj, row::Int)
#     new_obj = ConicObj()
#     for (var, coeff) in c
#         x1 = coeff[1][row, :]
#         x2 = coeff[2][row, :]
#         new_obj[var] = (x1,x2)
#     end
#     return new_obj
# end

# function Base.:*(v::Value, c::ConicObj)
#     # TODO: this part is time consuming, esp new_obj[var] = v * new_obj[var]...
#     new_obj = copy(c)
#     for var in keys(new_obj)
#         x1 = v * new_obj[var][1]
#         x2 = v * new_obj[var][2]
#         new_obj[var] = (x1,x2)
#     end
#     return new_obj
# end

# function promote_size(c::ConicObj, vectorized_size::Int)
#     new_obj = copy(c)
#     for var in keys(new_obj)
#         x1 = repeat(new_obj[var][1], vectorized_size, 1)
#         x2 = repeat(new_obj[var][2], vectorized_size, 1)
#         new_obj[var] = (x1,x2)
#     end
#     return new_obj
# end
