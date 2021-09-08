"""
    PositiveSemiDefinite(A::AbstractMatrix)

A wrapper used to indicate that a matrix is positive semi-definite
in order to possibly-expensive computations to verify that.

Currently only accepted by `quadform`.

!!! warning
    If the matrix is *not* positive semi-definite, you will likely
    get incorrect results.
"""
struct PositiveSemiDefinite{T<:AbstractMatrix}
    A::T
end
