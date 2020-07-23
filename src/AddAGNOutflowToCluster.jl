module AddAGNOutflowToCluster

    include("parameters.jl")
    include("io.jl")
    include("sample_particles.jl")

    export add_agn_outflow_to_halo,
           OutflowParameters

    function add_agn_outflow_to_halo(par::OutflowParameters)

        if par.verbose
            @info "Sampling positions"
            t1 = time_ns()
        end

        pos_outflow = sample_particles(par.Npart, par.theta_cone, par.r_max)

        if par.verbose
            t2 = time_ns()
            @info "Done!"
            @info "Took $(output_time(t1,t2)) s"
        end

        if par.verbose
            @info "Constructing outflow velocities"
            t1 = time_ns()
        end

        vel_outflow = construct_velocities(par.Npart, pos_outflow, par.v0)

        if par.verbose
            t2 = time_ns()
            @info "Done!"
            @info "Took $(output_time(t1,t2)) s"
        end

        if par.verbose
            @info "Constructing outflow B-field"
            t1 = time_ns()
        end

        B_outflow = construct_Bfield(par.Npart, pos_outflow, par.B0)

        if par.verbose
            t2 = time_ns()
            @info "Done!"
            @info "Took $(output_time(t1,t2)) s"
        end

        if par.verbose
            @info "Constructing outflow internal energy"
            t1 = time_ns()
        end

        u_outflow = construct_internal_energy(par)

        if par.verbose
            t2 = time_ns()
            @info "Done!"
            @info "Took $(output_time(t1,t2)) s"
        end

        if par.verbose
            t2 = time_ns()
            @info "Done!"
            @info "Took $(output_time(t1,t2)) s"
        end

        if par.verbose
            @info "Constructing outflow internal energy"
            t1 = time_ns()
        end

        m_outflow = construct_mass( par )

        if par.verbose
            t2 = time_ns()
            @info "Done!"
            @info "Took $(output_time(t1,t2)) s"
        end

        if par.verbose
            @info "Writing IC file"
            t1 = time_ns()
        end

        write_to_file( pos_outflow, vel_outflow, u_outflow, m_outflow, B_outflow, par)

        if par.verbose
            t2 = time_ns()
            @info "Done!"
            @info "Took $(output_time(t1,t2)) s"
        end

    end

end # module
