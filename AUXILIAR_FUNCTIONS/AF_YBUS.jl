# ===========================================
# Function to calculate the admittance matrix
# ===========================================
function Calculate_Ybus(DBUS::DataFrame, DCIR::DataFrame, nBUS::Int64, nCIR::Int64, base_MVA::Float64)
    # In MATPOWER and POWERMODELS, the TAP and SHIFT of the transformers are treated as "from" to "to"
    # This is the reason why they divide by TAP and not multiply when building the Ybus matrix

    Ybus = zeros(ComplexF64, nBUS, nBUS);   # Initialize the admittance matrix

    # Calculate the admittance matrix including the line data
    for i = 1:nCIR
        k = DCIR.from_bus[i] # Index from bus
        m = DCIR.to_bus[i]   # Index to bus
        if DCIR.l_status[i] == true

            ykm = 1 / (DCIR.l_res[i] + 1im*DCIR.l_reac[i]) # Series admittance
            bkm_sh = DCIR.l_sh_susp[i] / 2                 # Shunt admittance

            t_shift = deg2rad(DCIR.t_shift[i])

            Ybus[k,k] = Ybus[k,k] + ((1/DCIR.t_tap[i])^2 * ykm)+ (1im * bkm_sh)
            Ybus[k,m] = Ybus[k,m] - ((1/DCIR.t_tap[i]) * ykm * exp(1im * t_shift))
            Ybus[m,k] = Ybus[m,k] - ((1/DCIR.t_tap[i]) * ykm * exp(-1im * t_shift))
            Ybus[m,m] = Ybus[m,m] + ykm + (1im * bkm_sh)
        end
    end

    # Include the shunt components of the nodes
    for i = 1:nBUS
        Ybus[i,i] = Ybus[i,i] + (DBUS.g_sh[i] + 1im * DBUS.b_sh[i]) / base_MVA
    end

    return SparseArrays.sparse(Ybus)
end

# =========================================================================
# Function to calculate the admittance matrix for the "during fault" period
# =========================================================================
function Calculate_Ybus_fault(Ybus::SparseMatrixCSC, DBUS::DataFrame, DGEN::DataFrame, DGEN_DYN::DataFrame, nBUS::Int64, nGEN::Int64, base_MVA::Float64, bus_fault::Int64)

    # Considering the loads as constant admittances
    for bus in eachindex(DBUS.bus)
        Ybus[bus, bus] = Ybus[bus, bus] + ((DBUS.p_d[bus] - 1im*DBUS.q_d[bus]) / base_MVA)
    end

    # # Adding a high admittance in the faulted bus
    Ybus[bus_fault, bus_fault] = 1e6

    # Calculate susceptance of each generator
    nGEN_ON = count(isone, DGEN.g_status)
    indices_gen_ON = findall(!iszero, DGEN.g_status)

    Y_gen = zeros(ComplexF64, nGEN_ON)
    aux_count_gen_ON = 0
    for i in 1:nGEN
        if DGEN.g_status[i] == 1
            aux_count_gen_ON += 1
            Y_gen[aux_count_gen_ON] = 1 / (1im*DGEN_DYN.Xd_tr[i]) 
        end
    end

    # Construct incidence matrix for generators (like transmission lines from term_bus to gen_bus)
    A_gen = SparseArrays.sparse(
        vcat(DGEN_DYN.bus[indices_gen_ON], collect(1:nGEN_ON) .+ nBUS),   # row indices (terminal and internal buses)
        vcat(1:nGEN_ON, 1:nGEN_ON),                             # column indices
        vcat(ones(nGEN_ON), -ones(nGEN_ON)),                    # +1 at terminal bus, -1 at internal generator bus
        nBUS + nGEN_ON,                                         # total rows
        nGEN_ON                                                 # total columns
    )
    
    # Build the sparse diagonal matrix with the generators susceptances
    D = SparseArrays.spdiagm(Y_gen)

    # Compute B_gen_matrix = A_gen * diag(B_gen) * A_gen'
    Y_gen_matrix = A_gen * D * A_gen'

    # Create a larger sparse matrix to hold the expanded original B matrix
    Y_orig_expanded = SparseArrays.spzeros(ComplexF64, nBUS + nGEN_ON, nBUS + nGEN_ON)

    # Copy original matrix into the top-left corner
    Y_orig_expanded[1:nBUS, 1:nBUS] = Ybus

    # Now sum with the generator matrix
    Ybus_fault = Y_orig_expanded + Y_gen_matrix

    return SparseArrays.sparse(Ybus_fault)
end

# =======================================================================
# Function to calculate the admittance matrix for the "post-fault" period
# =======================================================================
function Calculate_Ybus_postf(DBUS::DataFrame, DCIR::DataFrame, DGEN::DataFrame, DGEN_DYN::DataFrame, nBUS::Int64, nCIR::Int64, nGEN::Int64, base_MVA::Float64, l_status::Vector{Int64})
    # In MATPOWER and POWERMODELS, the TAP and SHIFT of the transformers are treated as "from" to "to"
    # This is the reason why they divide by TAP and not multiply when building the Ybus matrix

    Ybus = zeros(ComplexF64, nBUS, nBUS);   # Initialize the admittance matrix

    # Calculate the admittance matrix including the line data
    for i = 1:nCIR
        k = DCIR.from_bus[i] # Index from bus
        m = DCIR.to_bus[i]   # Index to bus
        if l_status[i] == 1

            ykm = 1 / (DCIR.l_res[i] + 1im*DCIR.l_reac[i]) # Series admittance
            bkm_sh = DCIR.l_sh_susp[i] / 2                 # Shunt admittance

            t_shift = deg2rad(DCIR.t_shift[i])

            Ybus[k,k] = Ybus[k,k] + ((1/DCIR.t_tap[i])^2 * ykm)+ (1im * bkm_sh)
            Ybus[k,m] = Ybus[k,m] - ((1/DCIR.t_tap[i]) * ykm * exp(1im * t_shift))
            Ybus[m,k] = Ybus[m,k] - ((1/DCIR.t_tap[i]) * ykm * exp(-1im * t_shift))
            Ybus[m,m] = Ybus[m,m] + ykm + (1im * bkm_sh)
        end
    end

    # Include the shunt components of the nodes
    for bus = 1:nBUS
        Ybus[bus,bus] = Ybus[bus,bus] + (DBUS.g_sh[bus] + 1im * DBUS.b_sh[bus]) /base_MVA # Adding the shunt element

        Ybus[bus,bus] = Ybus[bus, bus] + ((DBUS.p_d[bus] - 1im*DBUS.q_d[bus]) / base_MVA)  # Adding the load as constant admittance
    end


    # Calculate susceptance of each generator
    nGEN_ON = count(isone, DGEN.g_status)
    indices_gen_ON = findall(!iszero, DGEN.g_status)

    Y_gen = zeros(ComplexF64, nGEN_ON)
    aux_count_gen_ON = 0
    for i in 1:nGEN
        if DGEN.g_status[i] == 1
            aux_count_gen_ON += 1
            Y_gen[aux_count_gen_ON] = 1 / (1im*DGEN_DYN.Xd_tr[i])
        end
    end

    # Construct incidence matrix for generators (like transmission lines from term_bus to gen_bus)
    A_gen = SparseArrays.sparse(
        vcat(DGEN_DYN.bus[indices_gen_ON], collect(1:nGEN_ON) .+ nBUS),   # row indices (terminal and internal buses)
        vcat(1:nGEN_ON, 1:nGEN_ON),                             # column indices
        vcat(ones(nGEN_ON), -ones(nGEN_ON)),                    # +1 at terminal bus, -1 at internal generator bus
        nBUS + nGEN_ON,                                         # total rows
        nGEN_ON                                                 # total columns
    )

    # Build the sparse diagonal matrix with the generators susceptances
    D = SparseArrays.spdiagm(Y_gen)

    # Compute B_gen_matrix = A_gen * diag(B_gen) * A_gen'
    Y_gen_matrix = A_gen * D * A_gen'

    # Create a larger sparse matrix to hold the expanded original B matrix
    Y_orig_expanded = SparseArrays.spzeros(ComplexF64, nBUS + nGEN_ON, nBUS + nGEN_ON)

    # Copy original matrix into the top-left corner
    Y_orig_expanded[1:nBUS, 1:nBUS] = SparseArrays.sparse(Ybus)

    # Now sum with the generator matrix
    Ybus_postf = Y_orig_expanded + Y_gen_matrix

    return SparseArrays.sparse(Ybus_postf)
end