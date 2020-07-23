using Base.Threads
using LinearAlgebra
using Statistics

function spherical_to_cartesian(r::Float64, ϕ::Float64, θ::Float64)

    x = r * cos(ϕ) * sin(θ)
    y = r * sin(ϕ) * sin(θ)
    z = r          * cos(θ)

    return x, y, z
end

function cartesian_to_spherical(x, y, z)
    
    r = sqrt(x^2 + y^2 + z^2)
    θ = atan(div(y, x))
    ϕ = acos(div(z,r))
    
    return r, θ, ϕ
end

function compute_outflow(pos::Vector{Float64}, prop0::Float64)

    # r, θ, ϕ = cartesian_to_spherical(pos[1], pos[2], pos[3])
    # return spherical_to_cartesian(prop0, θ, ϕ)
    
    #p = [0.0, 0.0, prop0] ⋅ pos 
    #return 0.0, 0.0, pos[3] * prop0
    return pos[1] * prop0, pos[2] * prop0, pos[3] * prop0
end

function sample_particles(Npart::Int64, θ_cone::Float64, r_max::Float64)

    θ_cone = deg2rad(θ_cone)

    Npart_half = floor(Int64, Npart*0.5)

    θ = zeros(Npart)
    θ[1:Npart_half] = 2.0 .* ( rand(Npart_half) .- 1.0 ) .* ( θ_cone * 0.5 )
    θ[Npart_half+1:end] = 2.0 .* ( rand(Npart_half) .- 1.0 ) .* ( θ_cone * 0.5 )
    ϕ = rand(Npart) .* 2π
    r = rand(Npart) .* r_max
    
    x = zeros(Npart)
    y = zeros(Npart)
    z = zeros(Npart)

    @threads for i = 1:Npart
        @inbounds x[i], y[i], z[i] = spherical_to_cartesian(r[i], ϕ[i], θ[i])
    end
    
    z[Npart_half+1:end] .*= -1.0

    return [x y z]
end


function construct_velocities(Npart::Int64, pos::Array{Float64,2}, v0::Float64)

    vx = zeros(Npart)
    vy = zeros(Npart)
    vz = zeros(Npart)

    for i = 1:Npart
        @inbounds vx[i], vy[i], vz[i] = compute_outflow(pos[i,:], v0)
    end

    return [vx vy vz]
end

function construct_Bfield(Npart::Int64, pos::Array{Float64,2}, B0::Float64)

    Bx = zeros(Npart)
    By = zeros(Npart)
    Bz = zeros(Npart)

    @threads for i = 1:Npart
        @inbounds Bx[i], By[i], Bz[i] = compute_outflow(pos[i,:], B0)
    end

    return [Bx By Bz]
end


function construct_internal_energy(par::OutflowParameters)
    return par.T .* ones(par.Npart) ./ par.units.T_K
end

function construct_mass(par::OutflowParameters)
    info_1d = Info_Line("POS", par.ic_format, 1, [1, 1, 1, 1, 1, 1])
    m    = read_block_by_name(par.input_snap, "MASS", info=info_1d, parttype=0)
    return ones(par.Npart) .* mean(m)
end