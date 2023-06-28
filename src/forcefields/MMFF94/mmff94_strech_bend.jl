@auto_hash_equals mutable struct MBend{T} <: AbstractForceFieldComponent{T}
    theta0   ::T
    thata_d  ::T
    theta    ::T
    ka       ::T
    at1      ::Atom{T}
    at2      ::Atom{T}
    at3      ::Atom{T}
    is_linear::Bool
    ATIJK    ::Int
    n1       ::Vector3{T}
    n2       ::Vector3{T}
end

@auto_hash_equals mutable struct MStretch{T} <: AbstractForceFieldComponent{T}
    at1      ::T
    at2      ::T
    kb       ::T
    t0       ::T
    delta_r  ::T
    sbmb     ::Bool
    empirical::Bool
    n        ::Vector3{T}
end

@auto_hash_equals mutable struct MStretchBend{T} <: AbstractForceFieldComponent{T}
    ff          ::ForceField{T}
    kba_ijk     ::T
    kba_kji     ::T
    sbtijk      ::T
    stretch_i_j ::MStretch
    stretch_j_k ::MStretch
    bend        ::MBend
end

function setup!(msb::MStretchBend{T}) where T<:Real
    
end