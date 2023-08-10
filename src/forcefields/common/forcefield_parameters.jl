export AbstractForceFieldParameters, extract_section

abstract type AbstractForceFieldParameters end

function extract_section(params::AbstractForceFieldParameters, section::String)
    if !haskey(params.sections, section)
        @error "extract_section(): Could not extract section $(section) from $(params.filename)!"
    end

    params.sections[section]
end

#adding unit mdyne = 1/1000 dyn
@unit mdyn    "mdyn"      Millidyne        (1/1000)*1Unitful.g*Unitful.cm/Unitful.s^2               true true

const force_prefactor = ustrip(u"kJ/mol/angstrom"/Constants.N_A |> u"N")
const joule_per_cal = ustrip(1u"cal" |> u"J")
const bend_k₀ = ustrip(0.021922u"cal" |> u"J")
const bend_k₀_forces = ustrip((bend_k₀)u"rad" |> u"°" )
const bend_kx = ustrip(143.9325u"cal" |> u"J")
const bend_k1 = -0.007
const stretch_k0 = ustrip((Float64(143.9325/2))u"cal" |> u"J")
const stretch_cubic_strength_constant = -2.
const stretch_kcs = 7. / 3.
const stretch_bend_k0 = ustrip(Float64(143.9325)u"cal" |> u"J") * 1u"°" |> u"rad" |> ustrip
const sb_constant = stretch_bend_k0/ (Constants.pi/180.0)  * force_prefactor
const k_torsion = 0.5u"cal" |> u"J" |> ustrip