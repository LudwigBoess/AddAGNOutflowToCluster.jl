using GadgetIO

function output_time(t1, t2)
    return Float64((t2-t1))*1.e-9
end

function read_data_gas(fi::String, dtype::DataType)

    info_3d = Info_Line("POS", dtype, 3, [1, 1, 1, 1, 1, 1])
    info_1d = Info_Line("POS", dtype, 1, [1, 1, 1, 1, 1, 1])

    pos  = read_block_by_name(fi, "POS",  info=info_3d, parttype=0)
    vel  = read_block_by_name(fi, "VEL",  info=info_3d, parttype=0) 
    u    = read_block_by_name(fi, "U",    info=info_1d, parttype=0)
    m    = read_block_by_name(fi, "MASS", info=info_1d, parttype=0)
    B    = read_block_by_name(fi, "BFLD", info=info_1d, parttype=0)

    return Float32.(pos), Float32.(vel), Float32.(u), Float32.(m), Float32.(B)
end

function read_data_collisionless(fi::String, parttype::Int64, dtype::DataType)

    info_3d = Info_Line("POS", dtype, 3, [1, 1, 1, 1, 1, 1])
    info_1d = Info_Line("POS", dtype, 1, [1, 1, 1, 1, 1, 1])

    pos  = read_block_by_name(fi, "POS",  info=info_3d, parttype=parttype)
    vel  = read_block_by_name(fi, "VEL",  info=info_3d, parttype=parttype)
    m    = read_block_by_name(fi, "MASS", info=info_1d, parttype=parttype) 

    return Float32.(pos), Float32.(vel), Float32.(m)
end

function correct_com(pos, par::OutflowParameters)

    x, v, m = read_data_collisionless(par.input_snap, 1, par.ic_format)

    com = zeros(3)

    @inbounds for i = 1:length(x[:,1])
        com[1] += x[i,1]
        com[2] += x[i,2]
        com[3] += x[i,3]
    end
    mtot = sum(m)

    com ./= mtot

    @inbounds for i = 1:length(pos[:,1])
        for j = 1:3
            pos[i,j] += com[j]
        end
    end

    return pos
end

function write_to_file(pos_outflow, vel_outflow, u_outflow, m_outflow, B_outflow, par::OutflowParameters)

    if par.verbose
        @info "  shifting outflow to COM"
        t1 = time_ns()
    end


    if par.verbose
        t2 = time_ns()
        @info "  Done!"
        @info "  Took $(output_time(t1,t2)) s"
    end

    if par.verbose
        @info "  Reading an merging snapshot and halo"
        t1 = time_ns()
    end

    # store number of particles in halo 
    Nhalo  = length(u_outflow)

    # read the header of the input snap
    header = head_to_obj(par.input_snap)

    # store total number of particles
    Ntotal = sum(header.npart) + Nhalo
    
    # set up id array
    ids = UInt32.(collect(1:Ntotal))

    # allocate empty arrays for particles
    pos = zeros(Float32, (Ntotal, 3) )
    vel = zeros(Float32, (Ntotal, 3) )
    u   = zeros(Float32, header.npart[1]+Nhalo )
    m   = zeros(Float32, Ntotal )
    B   = zeros(Float32, header.npart[1]+Nhalo )

    # assign halo particles
    pos[1:Nhalo,:] = Float32.(pos_outflow)
    vel[1:Nhalo,:] = Float32.(vel_outflow)
    u[1:Nhalo]     = Float32.(u_outflow)
    m[1:Nhalo]     = Float32.(m_outflow)
    B[1:Nhalo]     = Float32.(B_outflow)
 
    # read old gas particles
    Nstart = Nhalo + 1
    Npart  = header.npart[1]
    pos[Nstart:Nstart+Npart-1,:], vel[Nstart:Nstart+Npart-1,:], u[Nstart:Nstart+Npart-1], m[Nstart:Nstart+Npart-1], B[Nstart:Nstart+Npart-1] = read_data_gas(par.input_snap, par.ic_format)

    # read collisionless particles
    @inbounds for parttype = 1:5
        Npart = header.npart[parttype+1]
        if Npart > 0
            pos[Nstart:Nstart+Npart-1,:], vel[Nstart:Nstart+Npart-1,:], m[Nstart:Nstart+Npart-1] = read_data_collisionless(par.input_snap, parttype, par.ic_format)
            Nstart += Npart
        end
    end

    # update number of gas particles in snapshot
    header.npart[1] += Nhalo
    header.nall[1]  += Nhalo
    header.massarr  .= 0.0

    if par.verbose
        t2 = time_ns()
        @info "  Done!"
        @info "  Took $(output_time(t1,t2)) s"
    end

    if par.verbose
        @info "  Writing IC file"
        t1 = time_ns()
    end

    # write to file 
    f = open(par.output_file, "w")
    write_header(f, header)
    write_block(f, pos, "POS")
    write_block(f, vel, "VEL")
    write_block(f, ids, "ID")
    write_block(f, m, "MASS")
    write_block(f, u, "U")
    write_block(f, B, "BFLD")
    close(f)

    if par.verbose
        t2 = time_ns()
        @info "  Done!"
        @info "  Took $(output_time(t1,t2)) s"
    end

end