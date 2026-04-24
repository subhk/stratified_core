module PaperFiguresCommon

using LinearAlgebra

export YEAR,
       ModelParams,
       CouplingParams,
       kh,
       kz,
       k_total,
       omegaA,
       omegaI,
       omegaN,
       omegaM,
       gamma_mode,
       phi,
       uT_hat,
       raw_overlap,
       coupling,
       baseline_couplings,
       mode_spacing_ratio,
       rk4_envelope,
       diagnostics,
       exact_single_mode,
       peak_converted_fraction_exact,
       peak_converted_fraction_multimode

const YEAR = 365.25 * 24.0 * 3600.0

Base.@kwdef struct ModelParams
    Omega::Float64 = 7.3e-5
    vA::Float64 = 1.0e-2
    Ns::Float64 = 7.3e-5
    Hs::Float64 = 100.0e3
    eta::Float64 = 1.0
    lambda_h::Float64 = 2.0e6
end

Base.@kwdef struct CouplingParams
    lpen::Float64 = 50.0e3
    eps_overlap::Float64 = 6.0e-9
    gammaT::Float64 = 1.0 / (80.0 * YEAR)
    Nmodes::Int = 3
end

kh(p::ModelParams) = 2.0 * pi / p.lambda_h
kz(n::Integer, p::ModelParams) = (n + 0.5) * pi / p.Hs
k_total(n::Integer, p::ModelParams) = sqrt(kh(p)^2 + kz(n, p)^2)
omegaA(n::Integer, p::ModelParams) = p.vA * kz(n, p)
omegaI(n::Integer, p::ModelParams) = 2.0 * p.Omega * kz(n, p) / k_total(n, p)
omegaN(n::Integer, p::ModelParams) = p.Ns * kh(p) / k_total(n, p)
omegaM(n::Integer, p::ModelParams) = abs(omegaA(n, p)) * sqrt(omegaA(n, p)^2 + omegaN(n, p)^2) / abs(omegaI(n, p))
gamma_mode(n::Integer, p::ModelParams) = p.eta * kz(n, p)^2
phi(n::Integer, z::Real, p::ModelParams) = sqrt(2.0) * cos(kz(n, p) * (z + p.Hs))
uT_hat(z::Real, p::ModelParams, cp::CouplingParams) = exp(-(z + p.Hs) / cp.lpen)

function trapz(x::AbstractVector, y::AbstractVector)
    @assert length(x) == length(y)
    return sum((x[2:end] .- x[1:end-1]) .* (y[2:end] .+ y[1:end-1]) ./ 2.0)
end

function raw_overlap(n::Integer, p::ModelParams, cp::CouplingParams; nz::Int = 4001)
    z = collect(range(-p.Hs, 0.0, length = nz))
    integrand = [phi(n, zi, p) * uT_hat(zi, p, cp) for zi in z]
    return trapz(z, integrand) / p.Hs
end

function coupling(n::Integer, p::ModelParams, cp::CouplingParams, omegaT::Float64; nz::Int = 4001)
    return (p.Ns^2 / omegaT) * cp.eps_overlap * raw_overlap(n, p, cp; nz = nz)
end

function baseline_couplings(p::ModelParams, cp::CouplingParams, omegaT::Float64; nz::Int = 4001)
    return [coupling(n, p, cp, omegaT; nz = nz) for n in 0:(cp.Nmodes - 1)]
end

mode_spacing_ratio(n::Integer, p::ModelParams) = (omegaM(n + 1, p) - omegaM(n, p)) / gamma_mode(n, p)

function rhs_envelope(state::Vector{ComplexF64}, t::Float64,
                      gs::Vector{Float64}, deltas::Vector{Float64},
                      gammas::Vector{Float64}, gammaT::Float64)
    Nm = length(gs)
    dstate = zeros(ComplexF64, Nm + 1)
    AT = state[1]
    dstate[1] = -gammaT * AT + im * sum(gs[j] * state[j + 1] * exp(im * deltas[j] * t) for j in 1:Nm)
    for j in 1:Nm
        dstate[j + 1] = -gammas[j] * state[j + 1] + im * gs[j] * AT * exp(-im * deltas[j] * t)
    end
    return dstate
end

function rk4_envelope(p::ModelParams, cp::CouplingParams;
                      omegaT::Float64 = omegaM(0, p),
                      tmax::Float64 = 250.0 * YEAR,
                      dt::Float64 = 0.05 * YEAR,
                      A0::ComplexF64 = 1.0 + 0.0im,
                      nz::Int = 4001)
    gs = baseline_couplings(p, cp, omegaT; nz = nz)
    deltas = [omegaT - omegaM(n, p) for n in 0:(cp.Nmodes - 1)]
    gammas = [gamma_mode(n, p) for n in 0:(cp.Nmodes - 1)]

    times = collect(0.0:dt:tmax)
    A = zeros(ComplexF64, length(times), cp.Nmodes + 1)
    A[1, 1] = A0

    for i in 1:(length(times) - 1)
        y = vec(copy(A[i, :]))
        t = times[i]

        k1 = rhs_envelope(y, t, gs, deltas, gammas, cp.gammaT)
        k2 = rhs_envelope(y .+ 0.5 * dt .* k1, t + 0.5 * dt, gs, deltas, gammas, cp.gammaT)
        k3 = rhs_envelope(y .+ 0.5 * dt .* k2, t + 0.5 * dt, gs, deltas, gammas, cp.gammaT)
        k4 = rhs_envelope(y .+ dt .* k3, t + dt, gs, deltas, gammas, cp.gammaT)

        A[i + 1, :] .= y .+ dt .* (k1 .+ 2.0 .* k2 .+ 2.0 .* k3 .+ k4) ./ 6.0
    end

    return times, A, gs, deltas, gammas
end

function diagnostics(times::AbstractVector{<:Real}, A::AbstractMatrix{ComplexF64},
                     gammaT::Float64, gammas::Vector{Float64})
    ET = abs2.(A[:, 1])
    En = [abs2.(A[:, j]) for j in 2:size(A, 2)]

    Etot = copy(ET)
    for Ej in En
        Etot .+= Ej
    end

    sink = 2.0 .* gammaT .* ET
    for j in eachindex(En)
        sink .+= 2.0 .* gammas[j] .* En[j]
    end

    dissipation_integral = zeros(Float64, length(times))
    for i in 2:length(times)
        dt = times[i] - times[i - 1]
        dissipation_integral[i] = dissipation_integral[i - 1] + 0.5 * dt * (sink[i] + sink[i - 1])
    end

    epsE = (Etot .- Etot[1]) .+ dissipation_integral
    return ET, En, Etot, epsE
end

function exact_single_mode(times::AbstractVector{<:Real}, p::ModelParams, cp::CouplingParams;
                           delta::Float64 = 0.0,
                           A0::ComplexF64 = 1.0 + 0.0im,
                           nz::Int = 4001)
    omegaT = omegaM(0, p) + delta
    g0 = coupling(0, p, cp, omegaT; nz = nz)
    gamma0 = gamma_mode(0, p)
    chi = sqrt(complex(g0^2 - 0.25 * (cp.gammaT - gamma0 - im * delta)^2))

    prefactor = exp.(-0.5 .* (cp.gammaT + gamma0 + im * delta) .* times)
    AT = A0 .* prefactor .* (
        cos.(chi .* times) .+
        ((gamma0 - cp.gammaT + im * delta) / (2.0 * chi)) .* sin.(chi .* times)
    )
    B0 = im .* A0 .* (g0 / chi) .* prefactor .* sin.(chi .* times)
    return AT, B0, g0, gamma0
end

function peak_converted_fraction_exact(p::ModelParams, cp::CouplingParams;
                                       delta::Float64 = 0.0,
                                       tmax::Float64 = 250.0 * YEAR,
                                       dt::Float64 = 0.05 * YEAR,
                                       nz::Int = 4001)
    times = collect(0.0:dt:tmax)
    _, B0, _, _ = exact_single_mode(times, p, cp; delta = delta, nz = nz)
    return maximum(abs2.(B0))
end

function peak_converted_fraction_multimode(p::ModelParams, cp::CouplingParams;
                                           delta::Float64 = 0.0,
                                           tmax::Float64 = 250.0 * YEAR,
                                           dt::Float64 = 0.05 * YEAR,
                                           nz::Int = 4001)
    omegaT = omegaM(0, p) + delta
    _, A, _, _, _ = rk4_envelope(p, cp; omegaT = omegaT, tmax = tmax, dt = dt, nz = nz)
    return maximum(abs2.(A[:, 2]))
end

end # module
