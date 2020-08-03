using GadgetUnits

struct OutflowParameters

    Npart::Int64
    input_snap::String
    output_file::String
    v0::Float64
    B0::Float64
    T::Float64
    m_outflow::Float64
    center::Vector{Float64}
    theta_cone::Float64
    r_max::Float64
    verbose::Bool
    ic_format::DataType
    units::GadgetPhysical

    function OutflowParameters(Npart::Int64, 
                               input_snap::String, output_file::String, 
                               v0::Float64, B0::Float64,
                               T::Float64, m_outflow::Float64,
                               center::Array{Float64}=zeros(3),
                               theta_cone::Float64=5.0,
                               r_max::Float64=50.0
                               ;
                               verbose::Bool=true,
                               ic_format::DataType=Float32,
                               units::GadgetPhysical=GadgetPhysical()) 

        new(Npart, input_snap, output_file, v0, B0, T, m_outflow,
            center, theta_cone, r_max, verbose, ic_format,
            units)
    end
end