# ===================================================================================
#                  PRINT THE INPUT PARAMETERS IN TXT FILE
# ===================================================================================
function Print_Input_Parameters(case::String, 
    base_MVA::Float64, 
    δ_tol::Tuple{Float64, Float64},
    f_syn::Float64,
    bus_fault::Int,
    circ_trip::Vector,
    t_start_sim::Float64,
    t_end_sim::Float64,
    t_step::Float64,
    t_start_fault::Float64,
    clearing_time::Float64,
    t_clear_fault::Float64,
    current_path_folder::String,
    path_folder_results::String
    )

    cd(path_folder_results)

    filename = "input_parameters.txt"
    open(filename, "w") do io
        println(io, "******** Simulation Input Parameters ********")
        println(io, "==============================================")
        println(io, "Case:                         $case")
        println(io, "Base MVA:                     $base_MVA")
        println(io, "δ tolerance (deg):            [$(rad2deg(δ_tol[1])); +$(rad2deg(δ_tol[2]))]")
        println(io, "Synchronous frequency (Hz):   $f_syn")
        println(io, "Fault bus:                    $bus_fault")
        println(io, "Tripped circuits:             $(join(string.(circ_trip), ", "))")
        println(io, "Simulation start time (s):    $t_start_sim")
        println(io, "Simulation end time (s):      $t_end_sim")
        println(io, "Time step (s):                $t_step")
        println(io, "Fault start time (s):         $t_start_fault")
        println(io, "Fault lasts (s):              $clearing_time")
        println(io, "Fault is cleared at time (s): $t_clear_fault")
        println(io, "==============================================")
    end

    println("Input parameters successfully saved in: ", path_folder_results)

    cd(current_path_folder)
    
end

# ===================================================================================
#                   PRINT THE OPTIMIZATION MODEL IN TXT FILE
# ===================================================================================

# Function to write the AC-OPF model in a txt file
function Export_ACOPF_Model(model::Model, 
    V::OrderedDict{Int64, VariableRef}, 
    θ::OrderedDict{Int64, VariableRef}, 
    P_g::OrderedDict{Int64, VariableRef}, 
    Q_g::OrderedDict{Int64, VariableRef}, 
    P_ik::OrderedDict{Int64, VariableRef}, 
    Q_ik::OrderedDict{Int64, VariableRef}, 
    P_ki::OrderedDict{Int64, VariableRef}, 
    Q_ki::OrderedDict{Int64, VariableRef}, 
    eq_const_angle_sw::ConstraintRef, 
    eq_const_p_balance::OrderedDict{Int64, ConstraintRef}, 
    eq_const_q_balance::OrderedDict{Int64, ConstraintRef}, 
    eq_const_p_ik::OrderedDict{Int64, ConstraintRef}, 
    eq_const_q_ik::OrderedDict{Int64, ConstraintRef}, 
    eq_const_p_ki::OrderedDict{Int64, ConstraintRef}, 
    eq_const_q_ki::OrderedDict{Int64, ConstraintRef}, 
    ineq_const_s_ik::OrderedDict{Int64, ConstraintRef},
    ineq_const_s_ki::OrderedDict{Int64, ConstraintRef}, 
    ineq_const_diff_ang::OrderedDict{Int64, ConstraintRef},
    current_path_folder::String, 
    path_folder_results::String
    )

    cd(joinpath(path_folder_results,"ACOPF")) # Load the path folder for results

    # Open the file for writing
    open("model_summary.txt", "w") do io
        # Print the model to the file
        show(io, model)
    end

    # Desired key order to print the variables
    vector_dict_var = [V, θ, P_g, Q_g, P_ik, Q_ik, P_ki, Q_ki]

    open("ACOPF_model_details.txt", "w") do io

        # ------------------
        # Objective fuction
        # ------------------
        println(io, "=========")
        println(io, "Objective ")
        println(io, "=========")
        println(io, model.ext[:objective])
        println(io, "\n")

        # ---------------------------
        # Variables used in the model
        # ---------------------------
        println(io, "=========")
        println(io, "Variables")
        println(io, "=========")
        for i in eachindex(vector_dict_var)
            for (j, info) in vector_dict_var[i]
                println(io, "$j: ", info)
            end
        end
        println(io, "\n")

        # --------------------
        # Equality constraint
        # --------------------
        println(io, "===============================")
        println(io, "Equality Constraint Angle Swing ")
        println(io, "===============================")
        println(io, "1: ", eq_const_angle_sw)
        println(io, "\n")

        println(io, "===================================================")
        println(io, "Equality Constraints Active Power Balance for Buses ")
        println(io, "===================================================")
        for (i, info) in eq_const_p_balance 
            println(io, "$i: ", info) 
        end
        println(io, "\n")

        println(io, "=====================================================")
        println(io, "Equality Constraints Reactive Power Balance for Buses ")
        println(io, "=====================================================")
        for (i, info) in eq_const_q_balance 
            println(io, "$i: ", info) 
        end
        println(io, "\n")

        println(io, "=======================================================")
        println(io, "Equality Constraints Active Power Flow from Line i to k ")
        println(io, "=======================================================")
        for (i, info) in eq_const_p_ik 
            println(io, "$i: ", info)
        end
        println(io, "\n")

        println(io, "=========================================================")
        println(io, "Equality Constraints Reactive Power Flow from Line i to k ")
        println(io, "=========================================================")
        for (i, info) in eq_const_q_ik 
            println(io, "$i: ", info)
        end
        println(io, "\n")

        println(io, "=======================================================")
        println(io, "Equality Constraints Active Power Flow from Line k to i ")
        println(io, "=======================================================")
        for (i, info) in eq_const_p_ki 
            println(io, "$i: ", info)
        end
        println(io, "\n")

        println(io, "=========================================================")
        println(io, "Equality Constraints Reactive Power Flow from Line k to i ")
        println(io, "=========================================================")
        for (i, info) in eq_const_q_ki 
            println(io, "$i: ", info)
        end
        println(io, "\n")

        # ----------------------
        # Inequality constraint
        # ----------------------
        println(io, "===========================================================")
        println(io, "Inequality Constraints Capacity Power Flow from Line i to k ")
        println(io, "===========================================================")
        for (i, info) in ineq_const_s_ik
            println(io, "$i: ", info)
        end
        println(io, "\n")

        println(io, "===========================================================")
        println(io, "Inequality Constraints Capacity Power Flow from Line k to i ")
        println(io, "===========================================================")
        for (i, info) in ineq_const_s_ki
            println(io, "$i: ", info)
        end
        println(io, "\n")

        println(io, "==============================================================")
        println(io, "Inequality Constraints Voltage Angle Differences between Buses ")
        println(io, "==============================================================")
        for (i, info) in ineq_const_diff_ang
            println(io, "$i: ", info)
        end
        println(io, "\n")

        # Dicts with variables
        dicts_of_vars = [V, P_g, Q_g, P_ik, Q_ik, P_ki, Q_ki]  # OrderedDict{Int,VariableRef} or similar
        all_vars = Set(v for d in dicts_of_vars for v in values(d))

        println(io, "=========================================================")
        println(io, "Inequality Constraints Inferior Limits Decision Variables ")
        println(io, "=========================================================")
        const_lim_inf_decision_var = JuMP.all_constraints(model, VariableRef, MOI.LessThan{Float64})

        for (i, cref) in enumerate(const_lim_inf_decision_var)
            c_obj = JuMP.constraint_object(cref)  # returns ScalarConstraint
            var = c_obj.func                      # here .func is the VariableRef
            if var in all_vars
                println(io, "$i: ", cref)
            end
        end
        println(io, "\n")

        println(io, "=========================================================")
        println(io, "Inequality Constraints Superior Limits Decision Variables ")
        println(io, "=========================================================")
        const_lim_sup_decision_var = JuMP.all_constraints(model, VariableRef, MOI.GreaterThan{Float64})

        for (i, cref) in enumerate(const_lim_sup_decision_var)
            c_obj = JuMP.constraint_object(cref)  # returns ScalarConstraint
            var = c_obj.func                      # here .func is the VariableRef
            if var in all_vars
                println(io, "$i: ", cref)
            end
        end   
        println(io, "\n")

    end
    cd(current_path_folder)

    println("AC-OPF Model successfully saved as TXT file in: ", joinpath(path_folder_results,"ACOPF"))

end

# Function to write the AC-OPF model in a txt file
function Export_Dyn_Model(model::Model, 
    E::OrderedDict{Int, VariableRef}, 
    δ::OrderedDict{Int, VariableRef}, 
    Pm::OrderedDict{Int, VariableRef}, 
    eq_const_Pm::OrderedDict{Int, ConstraintRef},
    Pe_tf::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}, 
    δ_tf::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}, 
    Δω_tf::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}, 
    δCOI_tf::OrderedDict{Int, VariableRef},
    ΔωCOI_tf::OrderedDict{Int, VariableRef},
    eq_const_δCOI_tf::OrderedDict{Int, ConstraintRef},
    eq_const_ΔωCOI_tf::OrderedDict{Int, ConstraintRef},
    eq_const_Pe_tf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}},
    eq_const_δ_tf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}},
    eq_const_Δω_tf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}},
    ineq_const_δ_COI_tf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}},
    Pe_tpf::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}, 
    δ_tpf::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}, 
    Δω_tpf::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}, 
    δCOI_tpf::OrderedDict{Int64, VariableRef},
    ΔωCOI_tpf::OrderedDict{Int64, VariableRef},
    eq_const_δCOI_tpf::OrderedDict{Int, ConstraintRef},
    eq_const_ΔωCOI_tpf::OrderedDict{Int, ConstraintRef},
    eq_const_Pe_tpf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}},
    eq_const_δ_tpf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}},
    eq_const_Δω_tpf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}},
    ineq_const_δ_COI_tpf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}},
    current_path_folder::String, 
    path_folder_results::String
    )

    cd(joinpath(path_folder_results,"Transient_Stability")) # Load the path folder for results

    # Desired key order to print the variables
    vector_dict_var_pref      = [E, δ, Pm]
    vector_dict_var_fault     = [Pe_tf, δ_tf, Δω_tf]
    vector_dict_var_fault_COI = [δCOI_tf, ΔωCOI_tf]
    vector_dict_var_postf     = [Pe_tpf, δ_tpf, Δω_tpf]
    vector_dict_var_postf_COI = [δCOI_tpf, ΔωCOI_tpf]

    open("dynamic_model_details.txt", "w") do io

        # ---------------------------
        # Variables used in the model
        # ---------------------------
        println(io, "=================================")
        println(io, "Variables in the Pre-Fault Period")
        println(io, "=================================")
        for i in eachindex(vector_dict_var_pref)
            for (j, info) in vector_dict_var_pref[i]
                println(io, "$j: ", info)
            end
        end
        println(io, "\n")

        println(io, "=============================")
        println(io, "Variables in the Fault Period")
        println(io, "=============================")
        for i in eachindex(vector_dict_var_fault)
            for (j, infoj) in vector_dict_var_fault[i]
                for (k, infok) in infoj
                    println(io, "$k: ", infok)
                end
            end
        end
        for i in eachindex(vector_dict_var_fault_COI)
            for (j, info) in vector_dict_var_fault_COI[i]
                println(io, "$j: ", info)
            end
        end
        println(io, "\n")

        println(io, "==================================")
        println(io, "Variables in the Post-Fault Period")
        println(io, "==================================")
        for i in eachindex(vector_dict_var_postf)
            for (j, infoj) in vector_dict_var_postf[i]
                for (k, infok) in infoj
                    println(io, "$k: ", infok)
                end
            end
        end
        for i in eachindex(vector_dict_var_postf_COI)
            for (j, info) in vector_dict_var_postf_COI[i]
                println(io, "$j: ", info)
            end
        end
        println(io, "\n")

        # ---------------------
        # Equality constraints
        # ---------------------
        println(io, "====================================================")
        println(io, "Equality Constraints Mechanical Power of Generators ")
        println(io, "====================================================")
        for (i, info) in eq_const_Pm 
            println(io, "$i: ", info) 
        end
        println(io, "\n")


        println(io, "=====================================")
        println(io, "Equality Constraints Angle of the COI ")
        println(io, "=====================================")
        println(io, "------------")
        println(io, "Fault Period ")
        println(io, "------------")
        for (i, info) in eq_const_δCOI_tf 
            println(io, "$i: ", info) 
        end
        println(io, "------------------")
        println(io, "Post-Fault Period ")
        println(io, "------------------")
        for (i, info) in eq_const_δCOI_tpf 
            println(io, "$i: ", info) 
        end
        println(io, "\n")

        println(io, "================================================")
        println(io, "Equality Constraints Speed Deviation of the COI ")
        println(io, "================================================")
        println(io, "------------")
        println(io, "Fault Period ")
        println(io, "------------")
        for (i, info) in eq_const_ΔωCOI_tf 
            println(io, "$i: ", info) 
        end
        println(io, "-----------------")
        println(io, "Post-Fault Period ")
        println(io, "-----------------")
        for (i, info) in eq_const_ΔωCOI_tpf 
            println(io, "$i: ", info) 
        end
        println(io, "\n")

        println(io, "=====================================")
        println(io, "Equality Constraints Electrical Power ")
        println(io, "=====================================")
        println(io, "------------")
        println(io, "Fault Period ")
        println(io, "------------")
        for i in eachindex(eq_const_Pe_tf)
            for (j, info) in eq_const_Pe_tf[i]
                println(io, "$j: ", info)
            end
        end
        println(io, "------------------")
        println(io, "Post-Fault Period ")
        println(io, "------------------")
        for i in eachindex(eq_const_Pe_tpf)
            for (j, info) in eq_const_Pe_tpf[i]
                println(io, "$j: ", info)
            end
        end
        println(io, "\n")

        println(io, "==========================================")
        println(io, "Equality Constraints Angle Swing Equation")
        println(io, "==========================================")
        println(io, "------------")
        println(io, "Fault Period ")
        println(io, "------------")
        for i in eachindex(eq_const_δ_tf)
            for (j, info) in eq_const_δ_tf[i]
                println(io, "$j: ", info)
            end
        end
        println(io, "------------------")
        println(io, "Post-Fault Period ")
        println(io, "------------------")
        for i in eachindex(eq_const_δ_tpf)
            for (j, info) in eq_const_δ_tpf[i]
                println(io, "$j: ", info)
            end
        end
        println(io, "\n")

        println(io, "===================================================")
        println(io, "Equality Constraints Speed Deviation Swing Equation")
        println(io, "===================================================")
        println(io, "------------")
        println(io, "Fault Period ")
        println(io, "------------")
        for i in eachindex(eq_const_Δω_tf)
            for (j, info) in eq_const_Δω_tf[i]
                println(io, "$j: ", info)
            end
        end
        println(io, "------------------")
        println(io, "Post-Fault Period ")
        println(io, "------------------")
        for i in eachindex(eq_const_Δω_tpf)
            for (j, info) in eq_const_Δω_tpf[i]
                println(io, "$j: ", info)
            end
        end
        println(io, "\n")

        # ---------------------
        # Inequality constraints
        # ---------------------
        println(io, "===================================================")
        println(io, "Inequality Constraints Angle in Relation to the COI")
        println(io, "===================================================")
        println(io, "------------")
        println(io, "Fault Period ")
        println(io, "------------")
        for i in eachindex(ineq_const_δ_COI_tf)
            for (j, info) in ineq_const_δ_COI_tf[i]
                println(io, "$j: ", info)
            end
        end
        println(io, "------------------")
        println(io, "Post-Fault Period ")
        println(io, "------------------")
        for i in eachindex(ineq_const_δ_COI_tpf)
            for (j, info) in ineq_const_δ_COI_tpf[i]
                println(io, "$j: ", info)
            end
        end
        println(io, "\n")

        # Dicts with variables
        all_vars = Set(v for d in vector_dict_var_pref for v in values(d))

        println(io, "==============================================")
        println(io, "Constraints Superior Limits Decision Variables ")
        println(io, "==============================================")
        const_lim_inf_decision_var = JuMP.all_constraints(model, VariableRef, MOI.LessThan{Float64})

        for (i, cref) in enumerate(const_lim_inf_decision_var)
            c_obj = JuMP.constraint_object(cref)  # returns ScalarConstraint
            var = c_obj.func                      # here .func is the VariableRef
            if var in all_vars
                println(io, "$i: ", cref)
            end
        end
        println(io, "\n")

        println(io, "==============================================")
        println(io, "Constraints Inferior Limits Decision Variables ")
        println(io, "==============================================")
        const_lim_sup_decision_var = JuMP.all_constraints(model, VariableRef, MOI.GreaterThan{Float64})

        for (i, cref) in enumerate(const_lim_sup_decision_var)
            c_obj = JuMP.constraint_object(cref)  # returns ScalarConstraint
            var = c_obj.func                      # here .func is the VariableRef
            if var in all_vars
                println(io, "$i: ", cref)
            end
        end 
        println(io, "\n")  

    end
    cd(current_path_folder)

    println("Dynamic Model successfully saved as TXT file in: ", joinpath(path_folder_results,"Transient_Stability"))

end

# ===================================================================================
#                  PRINT THE REPORTS FOR THE AC OPF IN TXT AND CSV
# ===================================================================================

# Function to print and save variables according to the solution of the model
function Save_Solution_Model(model::Model, 
    V::OrderedDict{Int, VariableRef}, 
    θ::OrderedDict{Int, VariableRef}, 
    P_g::OrderedDict{Int, VariableRef}, 
    Q_g::OrderedDict{Int, VariableRef}, 
    P_ik::OrderedDict{Int, VariableRef}, 
    Q_ik::OrderedDict{Int, VariableRef}, 
    P_ki::OrderedDict{Int, VariableRef}, 
    Q_ki::OrderedDict{Int, VariableRef},
    bus_gen_circ_dict::OrderedDict,
    DBUS::DataFrame, 
    DGEN::DataFrame, 
    DGEN_DYN::DataFrame, 
    DCIR::DataFrame, 
    base_MVA::Float64, 
    nBUS::Int64, 
    nGEN::Int64, 
    nCIR::Int64, 
    bus_mapping::OrderedDict,
    reverse_bus_mapping::OrderedDict,
    current_path_folder::String,
    path_folder_results::String
    )


    println("===================================================")
    println("Objective Function: € "*string(round(JuMP.value.(model.ext[:objective]), digits=2))*"")
    println("===================================================")


    P_g_optim = [JuMP.value(v) for (i, v) in P_g]   # Get the results of the optimization process -> Variable P_g
    Q_g_optim = [JuMP.value(v) for (i, v) in Q_g]   # Get the results of the optimization process -> Variable Q_g
    S_g_optim = abs.(P_g_optim .+ 1im .* Q_g_optim) # Calculate the output complex power of the generator
    V_optim   = [JuMP.value(v) for (i, v) in V]     # Get the results of the optimization process -> Variable V
    θ_optim   = [JuMP.value(v) for (i, v) in θ]     # Get the results of the optimization process -> Variable θ

    Pik_optim = [JuMP.value(v) for (i, v) in P_ik] # Get the active power flow from bus i to bus k
    Qik_optim = [JuMP.value(v) for (i, v) in Q_ik] # Get the reactive power flow from bus i to bus k
    Pki_optim = [JuMP.value(v) for (i, v) in P_ki] # Get the active power flow from bus k to bus i
    Qki_optim = [JuMP.value(v) for (i, v) in Q_ki] # Get the reactive power flow from bus k to bus i

    Sik_optim = abs.(Pik_optim .+ 1im .* Qik_optim) # Apparent power flow from bus i to bus k
    Ski_optim = abs.(Pki_optim .+ 1im .* Qki_optim) # Apparent power flow from bus k to bus i

    # =============================================================================
    #                                   Generators
    # =============================================================================
    # Initialize Vector for all generators
    P_g_all = zeros(Float64, nGEN)
    Q_g_all = zeros(Float64, nGEN)
    S_g_all = zeros(Float64, nGEN)

    P_g_dict = Dict{Int, Float64}()
    Q_g_dict = Dict{Int, Float64}()
    S_g_dict = Dict{Int, Float64}()
    aux_count = 0
    for (i, id) in enumerate(DGEN.id)
        if DGEN.g_status[i] == 1
            aux_count += 1
            P_g_dict[id] = Float64(P_g_optim[aux_count])
            Q_g_dict[id] = Float64(Q_g_optim[aux_count])
            S_g_dict[id] = Float64(S_g_optim[aux_count])
        end
    end
    
    for (i, id) in enumerate(DGEN.id)
        if haskey(P_g_dict, id)
            P_g_all[i] = P_g_dict[id]  # Use optimized value
            Q_g_all[i] = Q_g_dict[id]  # Use optimized value
            S_g_all[i] = S_g_dict[id]
        end
    end

    # Calculate loading
    gen_loading_p = [DGEN.g_status[i] * ((P_g_all[i] - (DGEN.pg_min[i] / base_MVA)) / ((DGEN.pg_max[i] - DGEN.pg_min[i]) / base_MVA)) for i in 1:nGEN] # Generator loading -> Active power
    
    gen_loading_q = zeros(Float64, nGEN) # Generator loading -> Reactive power
    for i in 1:nGEN
        if Q_g_all[i] >= 0.0 # Check if the generator is providing reactive power to the system
            gen_loading_q[i] = DGEN.g_status[i] * (Q_g_all[i] / (DGEN.qg_max[i] / base_MVA))
        else # Or if the generator is consuming reactive power from the system
            gen_loading_q[i] = -DGEN.g_status[i] * (abs(Q_g_all[i]) / (abs(DGEN.qg_min[i]) / base_MVA))
        end
    end

    # =============================================================================
    #                                   Buses
    # =============================================================================
    # Defining a vector of Power Generated in each bus
    P_g_bus = zeros(Float64, nBUS) # Vector of active power generated at each bus
    Q_g_bus = zeros(Float64, nBUS) # Vector of reactive power generated at each bus
    for i in eachindex(DBUS.bus)
        indices_bus_gen = bus_gen_circ_dict[i][:gen_ids]
        if !isempty(indices_bus_gen)
            P_g_bus[i] = sum(P_g_all[indices_bus_gen])
            Q_g_bus[i] = sum(Q_g_all[indices_bus_gen])
        end
    end

    # =============================================================================
    #                                 Circuits
    # =============================================================================
    # Initialize Vector for all Circuits
    P_ik_all = zeros(Float64, nCIR)
    Q_ik_all = zeros(Float64, nCIR)
    S_ik_all = zeros(Float64, nCIR)
    P_ki_all = zeros(Float64, nCIR)
    Q_ki_all = zeros(Float64, nCIR)
    S_ki_all = zeros(Float64, nCIR)

    P_ik_dict = Dict{Int, Float64}()
    Q_ik_dict = Dict{Int, Float64}()
    S_ik_dict = Dict{Int, Float64}()
    P_ki_dict = Dict{Int, Float64}()
    Q_ki_dict = Dict{Int, Float64}()
    S_ki_dict = Dict{Int, Float64}()
    aux_count = 0
    for (i, id) in enumerate(DCIR.circ)
        if DCIR.l_status[i] == 1
            aux_count += 1
            P_ik_dict[id] = Float64(Pik_optim[aux_count])
            Q_ik_dict[id] = Float64(Qik_optim[aux_count])
            S_ik_dict[id] = Float64(Sik_optim[aux_count])
            P_ki_dict[id] = Float64(Pki_optim[aux_count])
            Q_ki_dict[id] = Float64(Qki_optim[aux_count])
            S_ki_dict[id] = Float64(Ski_optim[aux_count])
        end
    end

    for (i, id) in enumerate(DCIR.circ)
        if haskey(P_ik_dict, id)
            P_ik_all[i] = P_ik_dict[id]  # Use optimized value
            Q_ik_all[i] = Q_ik_dict[id]  # Use optimized value
            S_ik_all[i] = S_ik_dict[id]
            P_ki_all[i] = P_ki_dict[id]  # Use optimized value
            Q_ki_all[i] = Q_ki_dict[id]  # Use optimized value
            S_ki_all[i] = S_ki_dict[id]
        end
    end

    Plosses = P_ik_all + P_ki_all  # Active power losses in the branches
    Qlosses = Q_ik_all + Q_ki_all  # Reactive power losses in the branches
    

    # Calculate the loading of each branch according to the data provided in the input files
    # If there is no capacity limit defined in the input data file, the code assumes null loading
    all_cap = [DCIR.l_cap_1 DCIR.l_cap_2 DCIR.l_cap_3]
    circ_cap = zeros(Float64, nCIR)

    for i in 1:nCIR
        if any(!iszero, all_cap[i,:])
            index_cap = findfirst(!iszero, all_cap[i,:])
            circ_cap[i] = all_cap[i, index_cap]
        else
            circ_cap[i] = Inf
        end
    end
    circ_loading = [DCIR.l_status[lin] * abs(max(S_ik_all[lin], S_ki_all[lin])) / (circ_cap[lin] / base_MVA) for lin in eachindex(DCIR.circ)]

    # Correcting the buses labels
    bus_bus      = [reverse_bus_mapping[b] for b in DBUS.bus]
    gen_bus      = [reverse_bus_mapping[b] for b in DGEN.bus]
    gen_dyn_bus  = [reverse_bus_mapping[b] for b in DGEN_DYN.bus]
    from_bus     = [reverse_bus_mapping[b] for b in DCIR.from_bus]
    to_bus       = [reverse_bus_mapping[b] for b in DCIR.to_bus]

    # DataFrame to save the results related to the buses
    RBUS = DataFrame(
        bus  = bus_bus,                                             # Bus identifies  
        v    = V_optim,                                             # Voltage magnitude                            [p.u.]     
        θ    = round.(rad2deg.(θ_optim), digits=3),                 # Voltage angle                                [deg]  
        p    = round.((P_g_bus .* base_MVA) .- DBUS.p_d, digits=3), # Net Active power                             [MW] 
        q    = round.((Q_g_bus .* base_MVA) .- DBUS.q_d, digits=3), # Net Reactive power                           [MVAr] 
        p_g  = round.(P_g_bus .* base_MVA, digits=3),               # Active power generated                       [MW]  
        q_g  = round.(Q_g_bus .* base_MVA, digits=3),               # Reactive power generated                     [MVAr]  
        p_d  = round.(DBUS.p_d, digits=3),                          # Active power generated                       [MW]  
        q_d  = round.(DBUS.q_d, digits=3),                          # Reactive power demanded by load              [MVAr]  
        p_sh = round.(DBUS.g_sh .* (V_optim.^2), digits=3),         # Active power demanded by shunt conductance   [MW]   
        q_sh = round.(DBUS.b_sh .* (V_optim.^2), digits=3)          # Reactive power demanded by shunt suscpetance [MVAr]   
    )
    
    # DataFrame to save the results related to the circuits
    RCIR = DataFrame(
        circ      = DCIR.circ,                              # Circuit identifier
        from_bus  = from_bus,                               # From bus identifier
        to_bus    = to_bus,                                 # To bus identifier
        p_ik      = round.(P_ik_all .* base_MVA, digits=3), # Circuit active power flow from i to k   [MW]
        q_ik      = round.(Q_ik_all .* base_MVA, digits=3), # Circuit reactive power flow from i to k [MVAr]
        s_ik      = round.(S_ik_all .* base_MVA, digits=3), # Circuit apparent power flow from i to k [MVA]
        p_ki      = round.(P_ki_all .* base_MVA, digits=3), # Circuit active power flow from k to i   [MW]
        q_ki      = round.(Q_ki_all .* base_MVA, digits=3), # Circuit reactive power flow from k to i [MVAr]
        s_ki      = round.(S_ki_all .* base_MVA, digits=3), # Circuit apparent power flow from k to i [MVA]
        p_losses  = round.(Plosses .* base_MVA, digits=3),  # Losses of active power                  [MW]
        q_losses  = round.(Qlosses .* base_MVA, digits=3),  # Losses of reactive power                [MVAr]
        s_cap     = circ_cap,                               # Circuit maximum power capacity          [MVA]
        loading   = circ_loading                            # Circuit loading
    ) 

    # DataFrame to save the results related to the generators
    RGEN = DataFrame(
        id_gen    = DGEN.id,                                 # Generator ID
        id_bus    = gen_bus,                                 # Bus in which the generator is connected
        p_g       = round.((P_g_all .* base_MVA), digits=3), # Active power generated           [MW]
        q_g       = round.((Q_g_all .* base_MVA), digits=3), # Reactive power generated         [MVAr]
        s_g       = round.((S_g_all .* base_MVA), digits=3), # Apparent power generated         [MVA]
        loading_p = round.(gen_loading_p,         digits=3), # Generator active power loading
        loading_q = round.(gen_loading_q,         digits=3)  # Generator reactive power loading
    )

    Save_ResultsTXT_ACOPF(path_folder_results, RBUS, nBUS, RGEN, nGEN, RCIR, nCIR, Float64(JuMP.value.(model.ext[:objective]))) # Save the results in a TXT file
    Save_ResultsCSV_ACOPF(path_folder_results, RBUS, nBUS, RGEN, nGEN, RCIR, nCIR, Float64(JuMP.value.(model.ext[:objective]))) # Save the results in a CSV file

    cd(current_path_folder)

    return RBUS, RGEN, RCIR

end

# Save the reports of the power flow in TXT files
function Save_ResultsTXT_ACOPF(path_folder_results::String, RBUS::DataFrame, nBUS::Int64, RGEN::DataFrame, nGEN::Int64, RCIR::DataFrame, nCIR::Int64, objective::Float64)
    cd(joinpath(path_folder_results, "ACOPF"))

    # Summation of power generated and demmanded
    sum_Pg = 0.0
    sum_Qg = 0.0
    sum_Pd = 0.0
    sum_Qd = 0.0
    sum_Psh = 0.0
    sum_Qsh = 0.0
    sum_losses_P = 0.0
    sum_losses_Q = 0.0
    sum_gen_Pg = 0.0
    sum_gen_Qg = 0.0
    sum_gen_Sg = 0.0

    io = open("buses_report.txt", "w")
    @printf(io, "BUSES REPORT\n")
    @printf(io, "============================================================================================================================= \n")
    @printf(io, "   BUS    V (pu)     O (º)     P (MW)     Q (MVAr)      PG (MW)    QG (MVAr)    PD (MW)   QD (MVAr)     Psh (MW)   Qsh (MVAr) \n")
    @printf(io, "----------------------------------------------------------------------------------------------------------------------------- \n")
    for i = 1:nBUS
        @printf(io, " %4d    %6.4f    %6.2f    %8.2f    %8.2f    %8.2f    %8.2f    %8.2f    %8.2f    %8.2f    %8.2f\n", RBUS.bus[i], RBUS.v[i], RBUS.θ[i], RBUS.p[i], RBUS.q[i], RBUS.p_g[i], RBUS.q_g[i], RBUS.p_d[i], RBUS.q_d[i], RBUS.p_sh[i], RBUS.q_sh[i])
        sum_Pg += RBUS.p_g[i]
        sum_Qg += RBUS.q_g[i]
        sum_Pd += RBUS.p_d[i]
        sum_Qd += RBUS.q_d[i]
        sum_Psh += RBUS.p_sh[i]
        sum_Qsh += RBUS.q_sh[i]
    end
    @printf(io, "----------------------------------------------------------------------------------------------------------------------------- \n")
    @printf(io, " TOTAL:                                               %8.2f    %8.2f    %8.2f    %8.2f   %8.2f     %8.2f\n", sum_Pg, sum_Qg, sum_Pd , sum_Qd, sum_Psh , sum_Qsh)
    @printf(io, "============================================================================================================================= \n")
    @printf(io, "\n")
    close(io)  

    io = open("generators_report.txt", "w")
    @printf(io, "GENERATORS REPORT\n")
    @printf(io, "========================================================================== \n")
    @printf(io, "   ID     BUS     P_g (MW)  Q_g (MVAr)   S_g (MVA)   Loading_P   Loading_Q \n")
    @printf(io, "-------------------------------------------------------------------------- \n")
    for i = 1:nGEN
        @printf(io, " %4d   %4d    %8.2f    %8.2f    %8.2f    %8.4f    %8.4f   \n", RGEN.id_gen[i], RGEN.id_bus[i], RGEN.p_g[i], RGEN.q_g[i], RGEN.s_g[i], RGEN.loading_p[i], RGEN.loading_q[i])
        sum_gen_Pg += RGEN.p_g[i]
        sum_gen_Qg += RGEN.q_g[i]
        sum_gen_Sg += RGEN.s_g[i]
    end
    @printf(io, "-------------------------------------------------------------------------- \n")
    @printf(io, " TOTAL:         %8.4f   %8.4f   %8.4f\n", round(sum_gen_Pg, digits=3), round(sum_gen_Qg, digits=3), round(sum_gen_Sg, digits=3))
    @printf(io, "========================================================================== \n")
    close(io)
    
    io = open("circuits_report.txt", "w")
    @printf(io, "CIRCUITS REPORT\n")
    @printf(io, "=============================================================================================================================================== \n")
    @printf(io, "  CIRC    FROM    TO      Pik (MW)  Qik (MVAr)   Sik (MVA)   Pki (MW)  Qki (MVAr)    Ski (MVA)  Cap (MVA)    Loading    Ploss (MW)  Qloss (MVAr)\n")
    @printf(io, "----------------------------------------------------------------------------------------------------------------------------------------------- \n")
    for i = 1:nCIR
        @printf(io, " %4d   %4d   %4d    %8.2f    %8.2f    %8.2f    %8.2f   %8.2f     %8.2f   %8.4f     %8.4f     %8.3f     %8.3f\n", RCIR.circ[i], RCIR.from_bus[i], RCIR.to_bus[i], RCIR.p_ik[i], RCIR.q_ik[i], RCIR.s_ik[i], RCIR.p_ki[i], RCIR.q_ki[i], RCIR.s_ki[i], RCIR.s_cap[i], RCIR.loading[i], RCIR.p_losses[i], RCIR.q_losses[i])
        sum_losses_P += RCIR.p_losses[i]
        sum_losses_Q += RCIR.q_losses[i]
    end
    @printf(io, "----------------------------------------------------------------------------------------------------------------------------------------------- \n")
    @printf(io, " TOTAL:                                                                                                                  %8.3f     %8.3f\n", sum_losses_P , sum_losses_Q)
    @printf(io, "=============================================================================================================================================== \n")
    close(io)
    
    io = open("optimization_report.txt", "w")
    @printf(io, "OBJECTIVE\n")
    @printf(io, "================================ \n")
    @printf(io, "Total cost: (Euros) %8.2f \n", round(objective, digits = 2))
    @printf(io, "================================ \n")
    close(io)

    println("AC-OPF results successfully saved as TXT files in: ", joinpath(path_folder_results, "ACOPF"))
end

# Save the reports of the power flow in CSV files
function Save_ResultsCSV_ACOPF(path_folder_results::String, RBUS::DataFrame, nBUS::Int64, RGEN::DataFrame, nGEN::Int64, RCIR::DataFrame, nCIR::Int64, objective::Float64)
    cd(joinpath(path_folder_results, "ACOPF\\CSV"))
  
    # Save Bus Report as CSV
    df_buses = DataFrame(
        BUS = RBUS.bus,
        V_pu = round.(RBUS.v, digits = 2),
        Theta_deg = round.(RBUS.θ, digits = 2),
        P_MW = round.(RBUS.p, digits = 4),
        Q_MVAr = round.(RBUS.q, digits = 4),
        PG_MW = round.(RBUS.p_g, digits = 4),
        QG_MVAr = round.(RBUS.q_g, digits = 4),
        PD_MW = round.(RBUS.p_d, digits = 4),
        QD_MVAr = round.(RBUS.q_d, digits = 4),
        Psh_MW =  round.(RBUS.p_sh, digits = 4),
        Qsh_MVAr =  round.(RBUS.q_sh, digits = 4)
    )
    CSV.write("buses_report.csv", df_buses; delim=';', writeheader=true)

    # Save Generators Report as CSV
    df_generators = DataFrame(
        ID = RGEN.id_gen,
        BUS = RGEN.id_bus,
        P_MW = round.(RGEN.p_g, digits = 4),
        Q_MVAr = round.(RGEN.q_g, digits = 4),
        S_MVA = round.(RGEN.s_g, digits = 4),
        Loading_P = round.(RGEN.loading_p, digits = 4),
        Loading_Q = round.(RGEN.loading_q, digits = 4)
    )
    CSV.write("generators_report.csv", df_generators; delim=';', writeheader=true)

    # Save Circuit Report as CSV
    df_circuits = DataFrame(
        ID_CIRC = RCIR.circ,
        FROM_BUS = RCIR.from_bus,
        TO_BUS = RCIR.to_bus,
        Pik_MW = round.(RCIR.p_ik, digits = 4),
        Qik_MVAr = round.(RCIR.q_ik, digits = 4),
        Sik_MVA = round.(RCIR.s_ik, digits = 4),
        Pki_MW = round.(RCIR.p_ki, digits = 4),
        Qki_MVAr = round.(RCIR.q_ki, digits = 4),
        Ski_MVA = round.(RCIR.s_ki, digits = 4),
        Cap_MVA = round.(RCIR.s_cap, digits = 4),
        Loading = round.(RCIR.loading, digits = 4),
        Ploss_MW = round.(RCIR.p_losses, digits = 4),
        Qloss_MVAr = round.(RCIR.q_losses, digits = 4)
    )
    CSV.write("circuits_report.csv", df_circuits; delim=';', writeheader=true)

    # Save Optimization Report as CSV
    df_optimization = DataFrame(
        Metric = ["Total Cost (Euros)"],
        Value = [round(objective, digits=2)]
    )
    CSV.write("optimization_report.csv", df_optimization; delim=';', writeheader=true)

    println("AC-OPF results successfully saved as CSV files in: ", joinpath(path_folder_results, "ACOPF\\CSV"))
end

# ===================================================================================
#                  PRINT THE DUALS OF THE OPTIMIZATION PROBLEM
# ===================================================================================
# Function to obtain the duals of the AC-OPF
function Save_Duals_ACOPF_Model(model::Model,
    V::OrderedDict{Int, VariableRef}, 
    θ::OrderedDict{Int, VariableRef}, 
    P_g::OrderedDict{Int, VariableRef}, 
    Q_g::OrderedDict{Int, VariableRef}, 
    P_ik::OrderedDict{Int, VariableRef}, 
    Q_ik::OrderedDict{Int, VariableRef}, 
    P_ki::OrderedDict{Int, VariableRef}, 
    Q_ki::OrderedDict{Int, VariableRef},
    eq_const_angle_sw::ConstraintRef, 
    eq_const_p_balance::OrderedDict{Int64, ConstraintRef}, 
    eq_const_q_balance::OrderedDict{Int64, ConstraintRef}, 
    eq_const_p_ik::OrderedDict{Int64, ConstraintRef}, 
    eq_const_q_ik::OrderedDict{Int64, ConstraintRef}, 
    eq_const_p_ki::OrderedDict{Int64, ConstraintRef}, 
    eq_const_q_ki::OrderedDict{Int64, ConstraintRef}, 
    ineq_const_s_ik::OrderedDict{Int64, ConstraintRef},
    ineq_const_s_ki::OrderedDict{Int64, ConstraintRef}, 
    ineq_const_diff_ang::OrderedDict{Int64, ConstraintRef},
    base_MVA::Float64,
    current_path_folder::String,
    path_folder_results::String
    )

    cd(joinpath(path_folder_results, "ACOPF")) # Load the results path folder

    # -------------------------------------------------------------------------------------------------
    # Dual related to the equality constraint of angle at the swing bus
    dual_θ_SW = JuMP.dual.(eq_const_angle_sw) 

    # -------------------------------------------------------------------------------------------------
    # Dual related to the equality constraint active power balance
    dual_P_balance = [JuMP.dual(info) for (i, info) in eq_const_p_balance] ./ base_MVA

    # -------------------------------------------------------------------------------------------------
    # Dual related to the equality constraint reactive power balance
    dual_Q_balance = [JuMP.dual(info) for (i, info) in eq_const_q_balance] ./ base_MVA

    # -------------------------------------------------------------------------------------------------
    # Dual related to the equality constraint active power flow from i to k
    dual_Pik = [JuMP.dual(info) for (i, info) in eq_const_p_ik] ./ base_MVA

    # -------------------------------------------------------------------------------------------------
    # Dual related to the equality constraint reactive power flow from i to k
    dual_Qik = [JuMP.dual(info) for (i, info) in eq_const_q_ik] ./ base_MVA

    # -------------------------------------------------------------------------------------------------
    # Dual related to the equality constraint active power flow from k to i
    dual_Pki = [JuMP.dual(info) for (i, info) in eq_const_p_ki] ./ base_MVA

    # -------------------------------------------------------------------------------------------------
    # Dual related to the equality constraint reactive power flow from k to i
    dual_Qki = [JuMP.dual(info) for (i, info) in eq_const_q_ki] ./ base_MVA

    # -------------------------------------------------------------------------------------------------
    # Dual related to the inequality constraint apparent power flow capacity from i to k
    dual_Sik = [JuMP.dual(info) for (i, info) in ineq_const_s_ik] ./ base_MVA

    # -------------------------------------------------------------------------------------------------
    # Dual related to the inequality constraint apparent power flow capacity from k to i
    dual_Ski = [JuMP.dual(info) for (i, info) in ineq_const_s_ki] ./ base_MVA

    # -------------------------------------------------------------------------------------------------
    # Dual related to the inequality constraint voltage angle differences between buses
    dual_diff_ang = [JuMP.dual(info) for (i, info) in ineq_const_diff_ang]

    # -------------------------------------------------------------------------------------------------
    # Dual of the LOWER bound of the active power of each generator
    dual_LB_Pg = [JuMP.dual(LowerBoundRef(info)) for (i, info) in P_g] ./ base_MVA

    # Dual of the UPPER bound of the active power of each generator
    dual_UB_Pg = [JuMP.dual(UpperBoundRef(info)) for (i, info) in P_g] ./ base_MVA

    # -------------------------------------------------------------------------------------------------
    # Dual of the LOWER bound of the reactive power of each generator
    dual_LB_Qg = [JuMP.dual(LowerBoundRef(info)) for (i, info) in Q_g] ./ base_MVA

    # Dual of the UPPER bound of the reactive power of each generator
    dual_UB_Qg = [JuMP.dual(UpperBoundRef(info)) for (i, info) in Q_g] ./ base_MVA

    # -------------------------------------------------------------------------------------------------
    # Dual of the LOWER bound of each bus voltage
    dual_LB_V = [JuMP.dual(LowerBoundRef(info)) for (i, info) in V]

    # Dual of the UPPER bound of each bus voltage
    dual_UB_V = [JuMP.dual(UpperBoundRef(info)) for (i, info) in V]

    # -------------------------------------------------------------------------------------------------
    # Dual of the LOWER bound of the active power of circuits from -> to
    dual_LB_Pik = [JuMP.dual(LowerBoundRef(info)) for (i, info) in P_ik if JuMP.has_lower_bound(info)] ./ base_MVA

    # Dual of the UPPER bound of the active power of circuits from -> to
    dual_UB_Pik = [JuMP.dual(UpperBoundRef(info)) for (i, info) in P_ik if JuMP.has_upper_bound(info)] ./ base_MVA

    # Dual of the LOWER bound of the reactive power of circuits from -> to
    dual_LB_Qik = [JuMP.dual(LowerBoundRef(info)) for (i, info) in Q_ik if JuMP.has_lower_bound(info)] ./ base_MVA

    # Dual of the UPPER bound of the reactive power of circuits from -> to
    dual_UB_Qik = [JuMP.dual(UpperBoundRef(info)) for (i, info) in Q_ik if JuMP.has_upper_bound(info)] ./ base_MVA


    # -------------------------------------------------------------------------------------------------
    # Dual of the LOWER bound of the active power of circuits to -> from
    dual_LB_Pki = [JuMP.dual(LowerBoundRef(info)) for (i, info) in P_ki if JuMP.has_lower_bound(info)] ./ base_MVA

    # Dual of the UPPER bound of the active power of circuits to -> from
    dual_UB_Pki = [JuMP.dual(UpperBoundRef(info)) for (i, info) in P_ki if JuMP.has_upper_bound(info)] ./ base_MVA

    # Dual of the LOWER bound of the reactive power of circuits to -> from
    dual_LB_Qki = [JuMP.dual(LowerBoundRef(info)) for (i, info) in Q_ki if JuMP.has_lower_bound(info)] ./ base_MVA

    # Dual of the UPPER bound of the reactive power of circuits to -> from
    dual_UB_Qki = [JuMP.dual(UpperBoundRef(info)) for (i, info) in Q_ki if JuMP.has_upper_bound(info)] ./ base_MVA


    # ========== WRITE TO TXT FILE ==========
    open("ACOPF_duals.txt", "w") do io
        function write_dual_power(io, name, vec)
            println(io, "======================================")
            println(io, "          $name:")
            println(io, "======================================")
            for (i, val) in enumerate(vec)
                println(io, "[$i] =\t €/MW $val")
            end
            println(io)  # empty line between sections
        end

        function write_dual_others(io, name, vec)
            println(io, "======================================")
            println(io, "          $name:")
            println(io, "======================================")
            for (i, val) in enumerate(vec)
                println(io, "[$i] =\t $val")
            end
            println(io)  # empty line between sections
        end

        write_dual_others(io, "dual_θ_SW",     dual_θ_SW)
        write_dual_power(io, "dual_P_balance", dual_P_balance)
        write_dual_power(io, "dual_Q_balance", dual_Q_balance)
        write_dual_power(io, "dual_Pik",       dual_Pik)
        write_dual_power(io, "dual_Qik",       dual_Qik)
        write_dual_power(io, "dual_Pki",       dual_Pki)
        write_dual_power(io, "dual_Qki",       dual_Qki)
        write_dual_power(io, "dual_Sik",       dual_Sik)
        write_dual_power(io, "dual_Ski",       dual_Ski)
        write_dual_others(io, "dual_diff_ang", dual_diff_ang)

        write_dual_power(io,  "dual_LB_Pg", dual_LB_Pg)
        write_dual_power(io,  "dual_UB_Pg", dual_UB_Pg)
        write_dual_power(io,  "dual_LB_Qg", dual_LB_Qg)
        write_dual_power(io,  "dual_UB_Qg", dual_UB_Qg)
        write_dual_others(io, "dual_LB_V",  dual_LB_V)
        write_dual_others(io, "dual_UB_V",  dual_UB_V)

        write_dual_power(io, "dual_LB_Pik", dual_LB_Pik)
        write_dual_power(io, "dual_UB_Pik", dual_UB_Pik)
        write_dual_power(io, "dual_LB_Qik", dual_LB_Qik)
        write_dual_power(io, "dual_UB_Qik", dual_UB_Qik)

        write_dual_power(io, "dual_LB_Pki", dual_LB_Pki)
        write_dual_power(io, "dual_UB_Pki", dual_UB_Pki)
        write_dual_power(io, "dual_LB_Qki", dual_LB_Qki)
        write_dual_power(io, "dual_UB_Qki", dual_UB_Qki)
    end

    println("Duals of the AC-OPF model successfully saved as TXT file in: ", joinpath(path_folder_results, "ACOPF"))

    cd(current_path_folder)

end

# Function to obtain the duals of the Dynamic Simulation
function Save_Duals_Dynamic_Model(model::Model,
    E::OrderedDict{Int, VariableRef}, 
    δ::OrderedDict{Int, VariableRef}, 
    Pm::OrderedDict{Int, VariableRef}, 
    eq_const_Pm::OrderedDict{Int, ConstraintRef},
    Pe_tf::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}, 
    δ_tf::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}, 
    Δω_tf::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}, 
    δCOI_tf::OrderedDict{Int, VariableRef},
    ΔωCOI_tf::OrderedDict{Int, VariableRef},
    eq_const_δCOI_tf::OrderedDict{Int, ConstraintRef},
    eq_const_ΔωCOI_tf::OrderedDict{Int, ConstraintRef},
    eq_const_Pe_tf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}},
    eq_const_δ_tf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}},
    eq_const_Δω_tf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}},
    ineq_const_δ_COI_tf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}},
    Pe_tpf::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}, 
    δ_tpf::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}, 
    Δω_tpf::OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}, 
    δCOI_tpf::OrderedDict{Int64, VariableRef},
    ΔωCOI_tpf::OrderedDict{Int64, VariableRef},
    eq_const_δCOI_tpf::OrderedDict{Int, ConstraintRef},
    eq_const_ΔωCOI_tpf::OrderedDict{Int, ConstraintRef},
    eq_const_Pe_tpf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}},
    eq_const_δ_tpf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}},
    eq_const_Δω_tpf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}},
    ineq_const_δ_COI_tpf::OrderedDict{Int, OrderedDict{Int, ConstraintRef}},
    base_MVA::Float64,
    current_path_folder::String, 
    path_folder_results::String
    )

    cd(joinpath(path_folder_results, "Transient_Stability")) # Load the results path folder

    # -------------------------------------------------------------------------------------------------
    # Dual related to the equality constraint of mechanical power
    dual_Pm = [JuMP.dual(info) for (i, info) in eq_const_Pm] ./ base_MVA

    # -------------------------------------------------------------------------------------------------
    # Dual related to the equality constraint for the angle of the COI
    dual_δCOI_tf  = [JuMP.dual(info) for (i, info) in eq_const_δCOI_tf ] 
    dual_δCOI_tpf = [JuMP.dual(info) for (i, info) in eq_const_δCOI_tpf]
    dual_δCOI = vcat(dual_δCOI_tf, dual_δCOI_tpf)

    # -------------------------------------------------------------------------------------------------
    # Dual related to the equality constraint for the speed deviation of the COI
    dual_ΔωCOI_tf  = [JuMP.dual(info) for (i, info) in eq_const_ΔωCOI_tf ] 
    dual_ΔωCOI_tpf = [JuMP.dual(info) for (i, info) in eq_const_ΔωCOI_tpf]
    dual_ΔωCOI = vcat(dual_ΔωCOI_tf, dual_ΔωCOI_tpf)

    # -------------------------------------------------------------------------------------------------
    # Dual related to the equality constraint for the electrical power of each generator
    dual_Pe_tf  = OrderedDict(i => [JuMP.dual(v) for (_, v) in inner] ./ base_MVA for (i, inner) in eq_const_Pe_tf)
    dual_Pe_tpf = OrderedDict(i => [JuMP.dual(v) for (_, v) in inner] ./ base_MVA for (i, inner) in eq_const_Pe_tpf)

    dual_Pe = OrderedDict{Int, Vector{Float64}}()
    for k in keys(dual_Pe_tf)
        dual_Pe[k] = vcat(dual_Pe_tf[k], dual_Pe_tpf[k])
    end

    # -------------------------------------------------------------------------------------------------
    # Dual related to the equality constraint for the internal angle of each generator
    dual_δ_tf  = OrderedDict(i => [JuMP.dual(v) for (_, v) in inner] for (i, inner) in eq_const_δ_tf)
    dual_δ_tpf = OrderedDict(i => [JuMP.dual(v) for (_, v) in inner] for (i, inner) in eq_const_δ_tpf)

    dual_δ = OrderedDict{Int, Vector{Float64}}()
    for k in keys(dual_δ_tf)
        dual_δ[k] = vcat(dual_δ_tf[k], dual_δ_tpf[k])
    end

    # -------------------------------------------------------------------------------------------------
    # Dual related to the equality constraint for the speed deviation of each generator
    dual_Δω_tf  = OrderedDict(i => [JuMP.dual(v) for (_, v) in inner] for (i, inner) in eq_const_Δω_tf)
    dual_Δω_tpf = OrderedDict(i => [JuMP.dual(v) for (_, v) in inner] for (i, inner) in eq_const_Δω_tpf)

    dual_Δω = OrderedDict{Int, Vector{Float64}}()
    for k in keys(dual_Δω_tf)
        dual_Δω[k] = vcat(dual_Δω_tf[k], dual_Δω_tpf[k])
    end

    # -------------------------------------------------------------------------------------------------
    # Dual related to the inequality constraint for the internal angle of each generator in relation to the COI
    dual_δ_COI_tf  = OrderedDict(i => [JuMP.dual(v) for (_, v) in inner] for (i, inner) in ineq_const_δ_COI_tf)
    dual_δ_COI_tpf = OrderedDict(i => [JuMP.dual(v) for (_, v) in inner] for (i, inner) in ineq_const_δ_COI_tpf)

    dual_δ_COI = OrderedDict{Int, Vector{Float64}}()
    for k in keys(dual_δ_COI_tf)
        dual_δ_COI[k] = vcat(dual_δ_COI_tf[k], dual_δ_COI_tpf[k])
    end

    # -------------------------------------------------------------------------------------------------
    # Dual of the LOWER bound of the internal voltage of each generator
    dual_LB_E = [JuMP.dual(LowerBoundRef(info)) for (i, info) in E]

    # Dual of the UPPER bound of the internal voltage of each generator
    dual_UB_E = [JuMP.dual(UpperBoundRef(info)) for (i, info) in E]

    # -------------------------------------------------------------------------------------------------
    # Dual of the LOWER bound of the internal voltage of each generator
    dual_LB_δ = [JuMP.dual(LowerBoundRef(info)) for (i, info) in δ]

    # Dual of the UPPER bound of the internal voltage of each generator
    dual_UB_δ = [JuMP.dual(UpperBoundRef(info)) for (i, info) in δ]

    # -------------------------------------------------------------------------------------------------
    # Dual of the LOWER bound of the mechanical power of each generator
    dual_LB_Pm = [JuMP.dual(LowerBoundRef(info)) for (i, info) in Pm] ./ base_MVA

    # Dual of the UPPER bound of the mechanical power of each generator
    dual_UB_Pm = [JuMP.dual(UpperBoundRef(info)) for (i, info) in Pm] ./ base_MVA


    # ========== WRITE TO TXT FILE ==========
    open("dynamic_model_duals.txt", "w") do io
        function write_dual_power(io, name, vec)
            println(io, "======================================")
            println(io, "          $name:")
            println(io, "======================================")
            for (i, val) in enumerate(vec)
                println(io, "[$i] =\t €/MW $val")
            end
            println(io)  # empty line between sections
        end

        function write_dual_others(io, name, vec)
            println(io, "======================================")
            println(io, "          $name:")
            println(io, "======================================")
            for (i, val) in enumerate(vec)
                println(io, "[$i] =\t $val")
            end
            println(io)  # empty line between sections
        end

        write_dual_others(io, "dual_eq_Pm",     dual_Pm)
        write_dual_others(io, "dual_eq_δCOI",  dual_δCOI)
        write_dual_others(io, "dual_eq_ΔωCOI", dual_ΔωCOI)
        for (i, info) in dual_Pe
            write_dual_power(io, "dual_eq_Pe[G$i]", info)
        end

        for (i, info) in dual_δ
            write_dual_others(io, "dual_eq_δ[G$i]", info)
        end

        for (i, info) in dual_Δω
            write_dual_others(io, "dual_eq_Δω[G$i]", info)
        end

        dual_δ_COI
        for (i, info) in dual_δ_COI
            write_dual_others(io, "dual_ineq_δ_COI[G$i]", info)
        end

        write_dual_others(io, "dual_LB_E",  dual_LB_E)
        write_dual_others(io, "dual_UB_E",  dual_UB_E)
        write_dual_others(io, "dual_LB_δ",  dual_LB_δ)
        write_dual_others(io, "dual_UB_δ",  dual_UB_δ)
        write_dual_power(io,  "dual_LB_Pm", dual_LB_Pm)
        write_dual_power(io,  "dual_UB_Pm", dual_UB_Pm)
    end

    println("Duals of the dynamic model successfully saved as TXT file in: ", joinpath(path_folder_results, "Transient_Stability"))

    cd(current_path_folder)
    
end

# ===================================================================================
#                   MANAGE RESULTS FROM THE DYNAMIC SIMULATION
# ===================================================================================
# Function to save results across the whole time window (fault and post-fault)
function Manage_Dyn_Results(model::Model,
    DGEN_DYN::DataFrame,
    E::OrderedDict{Int, VariableRef}, 
    δ::OrderedDict{Int, VariableRef}, 
    Pm::OrderedDict{Int, VariableRef},
    Pe_tf::OrderedDict{Int, OrderedDict{Int, VariableRef}},
    δ_tf::OrderedDict{Int, OrderedDict{Int, VariableRef}},
    Δω_tf::OrderedDict{Int, OrderedDict{Int, VariableRef}},
    δCOI_tf::OrderedDict{Int, VariableRef},
    ΔωCOI_tf::OrderedDict{Int, VariableRef},
    Pe_tpf::OrderedDict{Int, OrderedDict{Int, VariableRef}},
    δ_tpf::OrderedDict{Int, OrderedDict{Int, VariableRef}},
    Δω_tpf::OrderedDict{Int, OrderedDict{Int, VariableRef}},
    δCOI_tpf::OrderedDict{Int, VariableRef},
    ΔωCOI_tpf::OrderedDict{Int, VariableRef},
    t_window_total::Vector{Float64},
    t_clear_fault::Float64,
    δ_tol::Tuple{Float64, Float64},
    base_MVA::Float64,
    f_syn::Float64,
    ω_syn::Float64,
    current_path_folder::String,
    path_folder_results::String
    )

    cd(joinpath(path_folder_results, "Transient_Stability"))

    E_values  = [JuMP.value(v) for (i, v) in E]
    δ_values  = [JuMP.value(v) for (i, v) in δ]
    Pm_values = [JuMP.value(v) for (i, v) in Pm] .* base_MVA

    # ====================================
    # Save short results in TXT
    # ====================================
    open("dynamic_model_short_results.txt", "w") do io
        # ---------------------------
        # Results for E, δ, and Pm
        # ---------------------------
        println(io, "==========================")
        println(io, "         E [p.u.]")
        println(io, "==========================")
        for (i, info) in E
            println(io, "$info: ", E_values[i])
        end
        println(io, "\n")

        println(io, "==========================")
        println(io, "         δ [deg]")
        println(io, "==========================")
        for (i, info) in δ
            println(io, "$info: ", rad2deg(δ_values[i]))
        end
        println(io, "\n")

        println(io, "==========================")
        println(io, "          Pm [MW]")
        println(io, "==========================")
        for (i, info) in δ
            println(io, "$info: ", Pm_values[i])
        end
        println(io, "\n")


    end

    # ====================================
    #                P_e
    # ====================================
    Pe_tf_values  = OrderedDict(i => [JuMP.value(v) for (_, v) in inner] .* base_MVA for (i, inner) in Pe_tf)  # Values of Pe during the fault
    Pe_tpf_values = OrderedDict(i => [JuMP.value(v) for (_, v) in inner] .* base_MVA for (i, inner) in Pe_tpf) # Values of Pe in the post-fault

    Pet = OrderedDict{Int, Vector{Float64}}()
    for k in keys(Pe_tf_values)
        Pet[k] = vcat(Pe_tf_values[k], Pe_tpf_values[k])
    end

    # ====================================
    #                 δ
    # ====================================
    δ_tf_values  = OrderedDict(i => [JuMP.value(v) for (_, v) in inner] for (i, inner) in δ_tf)   # Values of δ during the fault
    δ_tpf_values = OrderedDict(i => [JuMP.value(v) for (_, v) in inner] for (i, inner) in δ_tpf) # Values of δ in the post-fault

    δCOI_tf_values  = [JuMP.value(v) for (i, v) in δCOI_tf]  # Values of δ_COI during the fault
    δCOI_tpf_values = [JuMP.value(v) for (i, v) in δCOI_tpf] # Values of δ_COI in the post-fault

    δt = OrderedDict{Int, Vector{Float64}}()      # δ individual across the whole simulation
    δ_COIt = OrderedDict{Int, Vector{Float64}}()  # δ individual in relation to COI across the whole simulation

    δCOIt = vcat(δCOI_tf_values, δCOI_tpf_values) # δ do COI across the whole simulation

    for k in keys(δ_tf_values)
        δt[k]    = vcat(δ_tf_values[k], δ_tpf_values[k])
        δ_COIt[k] = δt[k] .- δCOIt
    end


    # ====================================
    #                 Δω
    # ====================================
    Δω_tf_values  = OrderedDict(i => [JuMP.value(v) for (_, v) in inner] for (i, inner) in Δω_tf)  # Values of Δω during the fault
    Δω_tpf_values = OrderedDict(i => [JuMP.value(v) for (_, v) in inner] for (i, inner) in Δω_tpf) # Values of Δω in the post-fault

    ΔωCOI_tf_values  = [JuMP.value(v) for (i, v) in ΔωCOI_tf]  # Values of Δω_COI during the fault
    ΔωCOI_tpf_values = [JuMP.value(v) for (i, v) in ΔωCOI_tpf] # Values of Δω_COI in the post-fault

    Δωt = OrderedDict{Int, Vector{Float64}}()     # Δω individual across the whole simulation
    Δω_COIt = OrderedDict{Int, Vector{Float64}}() # Δω individual in relation to COI across the whole simulation

    ΔωCOIt = vcat(ΔωCOI_tf_values, ΔωCOI_tpf_values) # Δω do COI across the whole simulation

    for k in keys(Δω_tf_values)
        Δωt[k]    = vcat(Δω_tf_values[k], Δω_tpf_values[k])
        Δω_COIt[k] = Δωt[k] .- ΔωCOIt
    end


    # ====================================
    #              RoCoF
    # ====================================
    time_RoCoF = []
    RoCoF = OrderedDict{Int, Vector{Float64}}()  # RoCoF individual across the whole simulation
    for (i, info) in Δωt
        time_RoCoF, RoCoF[i] = Calculate_Derivatives_CFD(t_window_total, f_syn .* (1.0 .+ info))
    end

    time_RoCoF, RoCoFCOI = Calculate_Derivatives_CFD(t_window_total, f_syn .* (1.0 .+ ΔωCOIt))

    # ====================================
    #          Kinetic Energy
    # ====================================
    Vke = OrderedDict{Int, Vector{Float64}}()  # Vke individual across the whole simulation
    for (i, info) in Δω_COIt
        Vke[i] = DGEN_DYN.H[i] .* ω_syn .* (info .^2)
    end

    # ====================================
    #        Accelerating Power
    # ====================================
    H_total   = sum(DGEN_DYN.H[Int.(keys(Pet))])

    Pacc      = OrderedDict{Int, Vector{Float64}}() # Pacc individual across the whole simulation
    PaccCOI   = OrderedDict{Int, Vector{Float64}}() # Pacc of the center of inertia across the whole simulation
    Pacc_COI  = OrderedDict{Int, Vector{Float64}}() # Pacc individual in relation to COI across the whole simulation
    aux_count = 0
    for (i, info) in Pet
        Pacc[i] = (Pm_values[i] .- info)

        aux_count += 1
        if aux_count == 1
            PaccCOI[1] = (Pm_values[i] .- info)
        else
            PaccCOI[1] = PaccCOI[1] .+ (Pm_values[i] .- info)
        end
    end
    for (i, info) in Pacc
        Pacc_COI[i] = (info .- ((DGEN_DYN.H[i] .* PaccCOI[1]) ./ H_total)) ./ base_MVA
    end
    
    # ====================================
    #          Potential Energy
    # ====================================
    Vpe = OrderedDict{Int, Vector{Float64}}() # Vpe individual across the whole simulation
    for k in keys(Pacc_COI)
        integral = []
        for (i, _) in enumerate(Pacc_COI[k])
            if i == 1
                push!(integral, 0.0)
            else
                push!(integral, -trapz(δ_COIt[k][1:i], Pacc_COI[k][1:i]))
            end
        end
        Vpe[k] = Float64.(integral) 
    end
    for (i, info) in Pacc_COI
        Pacc_COI[i] = info .* base_MVA
    end

    # ====================================
    #           SAVE FIGURES
    # ====================================
    Save_Dyn_Results_Plots(t_window_total, Pm_values, Pet, δt, δ_COIt, δCOIt, Δωt, Δω_COIt, ΔωCOIt, time_RoCoF, RoCoF, RoCoFCOI, Vke, Pacc, Pacc_COI, PaccCOI, Vpe, base_MVA, f_syn, δ_tol, t_clear_fault, current_path_folder, path_folder_results)

    # ====================================
    #           SAVE CSV FILES
    # ====================================
    
    Save_Dyn_Results_CSV(t_window_total, Pm_values, Pet, δt, δ_COIt, δCOIt, Δωt, Δω_COIt, ΔωCOIt, time_RoCoF, RoCoF, RoCoFCOI, Vke, Pacc, Pacc_COI, PaccCOI, Vpe, base_MVA, f_syn, current_path_folder, path_folder_results)

    cd(current_path_folder)

end

# Function to save the transient stability results in figures
function Save_Dyn_Results_Plots(t_window_total::Vector{Float64},
    Pm::Vector{Float64},
    Pet::OrderedDict{Int, Vector{Float64}},
    δt::OrderedDict{Int, Vector{Float64}},
    δ_COIt::OrderedDict{Int, Vector{Float64}},
    δCOIt::Vector{Float64},
    Δωt::OrderedDict{Int, Vector{Float64}},
    Δω_COIt::OrderedDict{Int, Vector{Float64}},
    ΔωCOIt::Vector{Float64},
    time_RoCoF::Vector{Float64},
    RoCoF::OrderedDict{Int, Vector{Float64}},
    RoCoFCOI::Vector{Float64},
    Vke::OrderedDict{Int, Vector{Float64}},
    Pacc::OrderedDict{Int, Vector{Float64}},
    Pacc_COI::OrderedDict{Int, Vector{Float64}},
    PaccCOI::OrderedDict{Int, Vector{Float64}},
    Vpe::OrderedDict{Int, Vector{Float64}},
    base_MVA::Float64,
    f_syn::Float64,
    δ_tol::Tuple{Float64, Float64},
    t_clear_fault::Float64,
    current_path_folder::String,
    path_folder_results::String
    )

    cd(joinpath(path_folder_results, "Transient_Stability"))

    index_t_ct = findfirst(x -> x >= t_clear_fault,  t_window_total) # Index of the vector when the time is equal or greater than the clearing time (only for plotting purposes)

    # ====================================
    #              Plots
    # ====================================

    #----------------------------------------
    # Plot the electrical power for each generator in MW
    plot_Pe = plot()  # create an empty plot

    for (gen_id, values) in Pet
        plot!(plot_Pe,
            t_window_total,
            values,
            lw = 3,
            label = "G$gen_id",
            ls=:solid
        )
    end
    plot_Pe = plot!(xlabel="t (s)", ylabel="P_e (MW)", title="Electrical Power",
    size=(1200,800), titlefont=font(40), xtickfont=font(35), ytickfont=font(35), guidefont=font(35), legendfont=font(15),
    fontfamily="Times New Roman", left_margin=10mm, bottom_margin=10mm, top_margin=10mm, right_margin=10mm,
    gridlinewidth=2, gridalpha=0.05, gridstyle=:dash)
    plot_Pe = plot!([t_clear_fault], seriestype=:vline, line=:dash, color=:magenta, label="t_ct")

    #----------------------------------------
    # Plot rotor angles in degrees
    plot_δ = plot()  # create an empty plot

    for (gen_id, values) in δt
        plot!(plot_δ,
            t_window_total,
            rad2deg.(values),
            lw = 3,
            label = "G$gen_id",
            ls=:solid
        )
    end
    plot_δ = plot!(t_window_total, rad2deg.(δCOIt)             , label="COI",   lw=2, ls=:solid, lc=:black)
    plot_δ = plot!(t_window_total, rad2deg.(δCOIt .+ δ_tol[1]) , label="- tol", lw=2, ls=:dash,  lc=:red)
    plot_δ = plot!(t_window_total, rad2deg.(δCOIt .+ δ_tol[2]) , label="+ tol", lw=2, ls=:dash,  lc=:red)
    plot_δ = plot!(xlabel="t (s)", ylabel="δ (deg)", title="Rotor Angles",
    size=(1200,800), titlefont=font(40), xtickfont=font(35), ytickfont=font(35), guidefont=font(35), legendfont=font(15),
    fontfamily="Times New Roman", left_margin=10mm, bottom_margin=10mm, top_margin=10mm, right_margin=10mm,
    gridlinewidth=2, gridalpha=0.05, gridstyle=:dash)
    plot_δ = plot!([t_clear_fault], seriestype=:vline, line=:dash, color=:magenta, label="t_ct")


    # Plot rotor angles in degrees in relation to the COI reference frame
    plot_δ_COI = plot()  # create an empty plot

    for (gen_id, values) in δ_COIt
        plot!(plot_δ_COI,
            t_window_total,
            rad2deg.(values),
            lw = 3,
            label = "G$gen_id",
            ls=:solid
        )
    end
    plot_δ_COI = plot!(xlabel="t (s)", ylabel="δ_COI (deg)", title="Rotor Angles (w.r.t. COI)",
    size=(1200,800), titlefont=font(40), xtickfont=font(35), ytickfont=font(35), guidefont=font(35), legendfont=font(15),
    fontfamily="Times New Roman", left_margin=10mm, bottom_margin=10mm, top_margin=10mm, right_margin=10mm,
    gridlinewidth=2, gridalpha=0.05, gridstyle=:dash)
    plot_δ_COI = plot!([t_clear_fault], seriestype=:vline, line=:dash, color=:magenta, label="t_ct")


    #----------------------------------------
    # Plot the speed deviation in p.u.
    plot_Δω = plot()  # create an empty plot

    for (gen_id, values) in Δωt
        plot!(plot_Δω,
            t_window_total,
            values,
            lw = 3,
            label = "G$gen_id",
            ls=:solid
        )
    end
    plot_Δω = plot!(t_window_total, ΔωCOIt              , label="COI",   lw=2, ls=:solid, lc=:black)
    plot_Δω = plot!(xlabel="t (s)", ylabel="Δω (p.u.)", title="Speed Deviation",
    size=(1200,800), titlefont=font(40), xtickfont=font(35), ytickfont=font(35), guidefont=font(35), legendfont=font(15),
    fontfamily="Times New Roman", left_margin=10mm, bottom_margin=10mm, top_margin=10mm, right_margin=10mm,
    gridlinewidth=2, gridalpha=0.05, gridstyle=:dash)
    plot_Δω = plot!([t_clear_fault], seriestype=:vline, line=:dash, color=:magenta, label="t_ct")


    # Plot speed deviation in relation to the COI reference frame
    plot_Δω_COI = plot()  # create an empty plot

    for (gen_id, values) in Δω_COIt
        plot!(plot_Δω_COI,
            t_window_total,
            values,
            lw = 3,
            label = "G$gen_id",
            ls=:solid
        )
    end
    plot_Δω_COI = plot!(xlabel="t (s)", ylabel="Δω_COI (p.u.)", title="Speed Deviation (w.r.t. COI)",
    size=(1200,800), titlefont=font(40), xtickfont=font(35), ytickfont=font(35), guidefont=font(35), legendfont=font(15),
    fontfamily="Times New Roman", left_margin=10mm, bottom_margin=10mm, top_margin=10mm, right_margin=10mm,
    gridlinewidth=2, gridalpha=0.05, gridstyle=:dash)
    plot_Δω_COI = plot!([t_clear_fault], seriestype=:vline, line=:dash, color=:magenta, label="t_ct")

    #--------------------------------------------
    # Plot the rotor speed in p.u.
    plot_Δω_1 = plot()  # create an empty plot

    for (gen_id, values) in Δωt
        plot!(plot_Δω_1,
            t_window_total,
            1 .+ values,
            lw = 3,
            label = "G$gen_id",
            ls=:solid
        )
    end
    plot_Δω_1 = plot!(t_window_total, 1 .+ (ΔωCOIt)              , label="COI",   lw=2, ls=:solid, lc=:black)
    plot_Δω_1 = plot!(xlabel="t (s)", ylabel="Δω (p.u.)", title="Rotor Speed",
    size=(1200,800), titlefont=font(40), xtickfont=font(35), ytickfont=font(35), guidefont=font(35), legendfont=font(15),
    fontfamily="Times New Roman", left_margin=10mm, bottom_margin=10mm, top_margin=10mm, right_margin=10mm,
    gridlinewidth=2, gridalpha=0.05, gridstyle=:dash)
    plot_Δω_1 = plot!([t_clear_fault], seriestype=:vline, line=:dash, color=:magenta, label="t_ct")


    # Plot rotor speed in p.u. in relation to the COI reference frame
    plot_Δω_COI_1 = plot()  # create an empty plot

    for (gen_id, values) in Δω_COIt
        plot!(plot_Δω_COI_1,
            t_window_total,
            1 .+ values,
            lw = 3,
            label = "G$gen_id",
            ls=:solid
        )
    end
    plot_Δω_COI_1 = plot!(xlabel="t (s)", ylabel="Δω_COI (p.u.)", title="Rotor Speed (w.r.t. COI)",
    size=(1200,800), titlefont=font(40), xtickfont=font(35), ytickfont=font(35), guidefont=font(35), legendfont=font(15),
    fontfamily="Times New Roman", left_margin=10mm, bottom_margin=10mm, top_margin=10mm, right_margin=10mm,
    gridlinewidth=2, gridalpha=0.05, gridstyle=:dash)
    plot_Δω_COI_1 = plot!([t_clear_fault], seriestype=:vline, line=:dash, color=:magenta, label="t_ct")

    #--------------------------------------------
    # Plot the frequency in Hz
    plot_f = plot()  # create an empty plot

    for (gen_id, values) in Δωt
        plot!(plot_f,
            t_window_total,
            f_syn .* (1 .+ values),
            lw = 3,
            label = "G$gen_id",
            ls=:solid
        )
    end

    plot_f = plot!(t_window_total, f_syn .* (1 .+ (ΔωCOIt))              , label="COI",   lw=2, ls=:solid, lc=:black)
    plot_f = plot!(xlabel="t (s)", ylabel="f (Hz)", title="Frequency",
    size=(1200,800), titlefont=font(40), xtickfont=font(35), ytickfont=font(35), guidefont=font(35), legendfont=font(15),
    fontfamily="Times New Roman", left_margin=10mm, bottom_margin=10mm, top_margin=10mm, right_margin=10mm,
    gridlinewidth=2, gridalpha=0.05, gridstyle=:dash)
    plot_f = plot!([t_clear_fault], seriestype=:vline, line=:dash, color=:magenta, label="t_ct")


    # Plot the frequency in relation to the COI in Hz
    plot_f_COI = plot()  # create an empty plot

    for (gen_id, values) in Δω_COIt
        plot!(plot_f_COI,
            t_window_total,
            f_syn .* (1 .+ values),
            lw = 3,
            label = "G$gen_id",
            ls=:solid
        )
    end
    plot_f_COI = plot!(xlabel="t (s)", ylabel="f (Hz)", title="Frequency (w.r.t. COI)",
    size=(1200,800), titlefont=font(40), xtickfont=font(35), ytickfont=font(35), guidefont=font(35), legendfont=font(15),
    fontfamily="Times New Roman", left_margin=10mm, bottom_margin=10mm, top_margin=10mm, right_margin=10mm,
    gridlinewidth=2, gridalpha=0.05, gridstyle=:dash)
    plot_f_COI = plot!([t_clear_fault], seriestype=:vline, line=:dash, color=:magenta, label="t_ct")

    #--------------------------------------------
    # Plot the RoCoF in Hz/s
    plot_RoCoF = plot()  # create an empty plot

    for (gen_id, values) in RoCoF
        plot!(plot_RoCoF,
            time_RoCoF,
            values,
            lw = 3,
            label = "G$gen_id",
            ls=:solid
        )
    end

    plot_RoCoF = plot!(time_RoCoF, RoCoFCOI, label="COI",   lw=2, ls=:solid, lc=:black)
    plot_RoCoF = plot!(xlabel="t (s)", ylabel="df/dt (Hz/sec)", title="RoCoF",
    size=(1200,800), titlefont=font(40), xtickfont=font(35), ytickfont=font(35), guidefont=font(35), legendfont=font(15),
    fontfamily="Times New Roman", left_margin=10mm, bottom_margin=10mm, top_margin=10mm, right_margin=10mm,
    gridlinewidth=2, gridalpha=0.05, gridstyle=:dash)
    plot_RoCoF = plot!([t_clear_fault], seriestype=:vline, line=:dash, color=:magenta, label="t_ct")


    #--------------------------------------------
    # Plot δ vs ω
    plot_δω = plot()  # create an empty plot

    for (gen_id, values) in δ_COIt
        plot!(plot_δω,
            rad2deg.(values),
            1.0 .+ Δω_COIt[gen_id],
            lw = 3,
            label = "G$gen_id",
            ls=:solid
        )
    end
    plot_δω = plot!(xlabel="δ (deg)", ylabel="Δω (p.u.)", title="Rotor speed vs Angle",
    size=(1200,800), titlefont=font(40), xtickfont=font(35), ytickfont=font(35), guidefont=font(35), legendfont=font(15),
    fontfamily="Times New Roman", left_margin=10mm, bottom_margin=10mm, top_margin=10mm, right_margin=10mm,
    gridlinewidth=2, gridalpha=0.05, gridstyle=:dash)

    for (gen_id, values) in δ_COIt
        scatter!(plot_δω,
            [rad2deg.(values[index_t_ct])],
            [1.0 .+ Δω_COIt[gen_id][index_t_ct]],
            shape=:circle, 
            msize= 7, 
            color=:black, 
            label=""
        )
    end

    #--------------------------------------------
    # Plot δ vs Pe and Pm
    plot_δ_COIPePm = plot()  # create an empty plot
 
    for (gen_id, values) in δ_COIt
        # Pe
        plot!(plot_δ_COIPePm,
            rad2deg.(values),
            Pet[gen_id],
            lw = 3,
            label = "Pe G$gen_id",
            ls=:solid
        )
        # Pm
        plot!(plot_δ_COIPePm,
            rad2deg.(values),
            Pm[gen_id] .* ones(Float64, length(values)),
            lw = 3,
            label = "Pm G$gen_id",
            ls=:dash
        )
    end
    
    plot_δ_COIPePm = plot!(xlabel="δ (deg)", ylabel="P (MW)", title="P vs Angle",
    size=(1200,800), titlefont=font(40), xtickfont=font(35), ytickfont=font(35), guidefont=font(35), legendfont=font(15),
    fontfamily="Times New Roman", left_margin=10mm, bottom_margin=10mm, top_margin=10mm, right_margin=10mm,
    gridlinewidth=2, gridalpha=0.05, gridstyle=:dash)

    #--------------------------------------------
    # Plot t vs Pacc 
    plot_Pacc = plot()  # create an empty plot
    # Pacc
    for (gen_id, values) in Pacc
        plot!(plot_Pacc,
            t_window_total,
            values,
            lw = 3,
            label = "G$gen_id",
            ls=:solid
        )
    end
    plot_Pacc = plot!(t_window_total, PaccCOI[1], label="COI", lw=2, ls=:solid, lc=:black)
    plot_Pacc = plot!(xlabel="t (s)", ylabel="P_acc (MW)", title="Accelerating Power",
    size=(1200,800), titlefont=font(40), xtickfont=font(35), ytickfont=font(35), guidefont=font(35), legendfont=font(15),
    fontfamily="Times New Roman", left_margin=10mm, bottom_margin=10mm, top_margin=10mm, right_margin=10mm,
    gridlinewidth=2, gridalpha=0.05, gridstyle=:dash)
    plot_Pacc = plot!([t_clear_fault], seriestype=:vline, line=:dash, color=:magenta, label="t_ct")

    #--------------------------------------------
    # Plot t vs Pacc (in relation to COI)
    plot_Pacc_COI = plot()  # create an empty plot
    # Pacc
    for (gen_id, values) in Pacc_COI
        plot!(plot_Pacc_COI,
            t_window_total,
            values,
            lw = 3,
            label = "G$gen_id",
            ls=:solid
        )
    end
    plot_Pacc_COI = plot!(xlabel="t (s)", ylabel="P_acc_COI (MW)", title="Accelerating Power (w.r.t. COI)",
    size=(1200,800), titlefont=font(40), xtickfont=font(35), ytickfont=font(35), guidefont=font(35), legendfont=font(15),
    fontfamily="Times New Roman", left_margin=10mm, bottom_margin=10mm, top_margin=10mm, right_margin=10mm,
    gridlinewidth=2, gridalpha=0.05, gridstyle=:dash)
    plot_Pacc_COI = plot!([t_clear_fault], seriestype=:vline, line=:dash, color=:magenta, label="t_ct")

    #--------------------------------------------
    # Plot t vs Vke
    plot_Vke = plot()  # create an empty plot
    # Vke
    for (gen_id, values) in Vke
        plot!(plot_Vke,
            t_window_total,
            values,
            lw = 3,
            label = "G$gen_id",
            ls=:solid
        )
    end
    plot_Vke = plot!(xlabel="t (s)", ylabel="V_ke (p.u. ⋅ rad)", title="Vke vs Time",
    size=(1200,800), titlefont=font(40), xtickfont=font(35), ytickfont=font(35), guidefont=font(35), legendfont=font(15),
    fontfamily="Times New Roman", left_margin=10mm, bottom_margin=10mm, top_margin=10mm, right_margin=10mm,
    gridlinewidth=2, gridalpha=0.05, gridstyle=:dash)
    plot_Vke = plot!([t_clear_fault], seriestype=:vline, line=:dash, color=:magenta, label="t_ct")


    #--------------------------------------------
    # Plot t vs Vpe
    plot_Vpe = plot()  # create an empty plot
    # Vpe
    for (gen_id, values) in Vpe
        plot!(plot_Vpe,
            t_window_total,
            values,
            lw = 3,
            label = "G$gen_id",
            ls=:solid
        )
    end
    plot_Vpe = plot!(xlabel="t (s)", ylabel="V_pe (p.u. ⋅ rad)", title="Vpe vs Time",
    size=(1200,800), titlefont=font(40), xtickfont=font(35), ytickfont=font(35), guidefont=font(35), legendfont=font(15),
    fontfamily="Times New Roman", left_margin=10mm, bottom_margin=10mm, top_margin=10mm, right_margin=10mm,
    gridlinewidth=2, gridalpha=0.05, gridstyle=:dash)
    plot_Vpe = plot!([t_clear_fault], seriestype=:vline, line=:dash, color=:magenta, label="t_ct")

    #--------------------------------------------
    # Plot t vs Ve
    plot_Ve = plot()  # create an empty plot
    # Ve
    for (gen_id, values) in Vpe
        plot!(plot_Ve,
            t_window_total,
            Vke[gen_id] .+ values,
            lw = 3,
            label = "G$gen_id",
            ls=:solid
        )
    end
    plot_Ve = plot!(xlabel="t (s)", ylabel="V_e (p.u. ⋅ rad)", title="Total Energy vs Time",
    size=(1200,800), titlefont=font(40), xtickfont=font(35), ytickfont=font(35), guidefont=font(35), legendfont=font(15),
    fontfamily="Times New Roman", left_margin=10mm, bottom_margin=10mm, top_margin=10mm, right_margin=10mm,
    gridlinewidth=2, gridalpha=0.05, gridstyle=:dash)
    plot_Ve = plot!([t_clear_fault], seriestype=:vline, line=:dash, color=:magenta, label="t_ct")

    # ====================================
    #            Save Plots
    # ====================================
    cd(joinpath(path_folder_results, "Transient_Stability\\Figures"))
    savefig(plot_Pe,        "Electrical Power vs time.svg")
    savefig(plot_δ,         "Delta vs time.svg")
    savefig(plot_δ_COI,     "Delta (ref COI) vs time.svg")
    savefig(plot_Δω,        "Speed Deviation vs time.svg")
    savefig(plot_Δω_COI,    "Speed Deviation (ref COI) vs time.svg")
    savefig(plot_Δω_1,      "Rotor speed vs time.svg")
    savefig(plot_Δω_COI_1,  "Rotor speed (ref COI) vs time.svg")
    savefig(plot_f,         "Frequency vs time.svg")
    savefig(plot_f_COI,     "Frequency (ref COI) vs time.svg")
    savefig(plot_RoCoF,     "RoCoF vs time.svg")
    savefig(plot_δω,        "Rotor speed vs delta.svg")
    savefig(plot_δ_COIPePm, "Power vs delta.svg")
    savefig(plot_Pacc,      "Accelerating power vs time.svg")
    savefig(plot_Pacc_COI,  "Accelerating power (ref COI) vs time.svg")
    savefig(plot_Vke,       "Kinetic Energy vs time.svg")
    savefig(plot_Vpe,       "Potential Energy vs time.svg")
    savefig(plot_Ve,        "Total Energy vs time.svg")

    cd(current_path_folder)

    println("Figures of the dynamic model successfully saved as PNG files in: ", joinpath(path_folder_results, "Transient_Stability\\Figures"))
end

# Function to save the transient stability results in CSV files
function Save_Dyn_Results_CSV(t_window_total::Vector{Float64},
    Pm::Vector{Float64},
    Pet::OrderedDict{Int, Vector{Float64}},
    δt::OrderedDict{Int, Vector{Float64}},
    δ_COIt::OrderedDict{Int, Vector{Float64}},
    δCOIt::Vector{Float64},
    Δωt::OrderedDict{Int, Vector{Float64}},
    Δω_COIt::OrderedDict{Int, Vector{Float64}},
    ΔωCOIt::Vector{Float64},
    time_RoCoF::Vector{Float64},
    RoCoF::OrderedDict{Int, Vector{Float64}},
    RoCoFCOI::Vector{Float64},
    Vke::OrderedDict{Int, Vector{Float64}},
    Pacc::OrderedDict{Int, Vector{Float64}},
    Pacc_COI::OrderedDict{Int, Vector{Float64}},
    PaccCOI::OrderedDict{Int, Vector{Float64}},
    Vpe::OrderedDict{Int, Vector{Float64}},
    base_MVA::Float64,
    f_syn::Float64,
    current_path_folder::String,
    path_folder_results::String
    )

    cd(joinpath(path_folder_results, "Transient_Stability"))
    
    # ========================================================
    # Names of the variables
    # ========================================================

    gen_names = ["G$(i)" for (i, inner) in Pet]

    delta_names = ["G$(i)" for (i, inner) in δt]
    push!(delta_names, "COI")

    omega_names = ["G$(i)" for (i, inner) in Δωt]
    push!(omega_names, "COI")

    # ========================================================
    # Creating the matrices
    # ========================================================
    # Electrical Power
    Pe_matrix = zeros(Float64, length(t_window_total), length(Pet))
    aux_count = 0
    for (gen_id, values) in Pet
        aux_count += 1
        Pe_matrix[:,aux_count] = values
    end
    df_Pe = DataFrame(hcat(t_window_total, Pe_matrix), vcat("t", gen_names))

    # Accelerating Power
    Pacc_matrix = zeros(Float64, length(t_window_total), length(Pacc)+1)
    aux_count = 0
    for (gen_id, values) in Pacc
        aux_count += 1
        Pacc_matrix[:,aux_count] = values
    end
    Pacc_matrix[:,end] = PaccCOI[1]
    df_Pacc = DataFrame(hcat(t_window_total, Pacc_matrix), vcat("t", delta_names))

    # Accelerating Power in relation to COI
    Pacc_COI_matrix = zeros(Float64, length(t_window_total), length(Pacc_COI))
    aux_count = 0
    for (gen_id, values) in Pacc_COI
        aux_count += 1
        Pacc_COI_matrix[:,aux_count] = values
    end
    df_Pacc_COI = DataFrame(hcat(t_window_total, Pacc_COI_matrix), vcat("t", gen_names))

    # Absolute values of δ
    δ_matrix = zeros(Float64, length(t_window_total), length(δt)+1)
    aux_count = 0
    for (gen_id, values) in δt
        aux_count += 1
        δ_matrix[:,aux_count] = rad2deg.(values)
    end
    δ_matrix[:,end] = rad2deg.(δCOIt)
    df_δ = DataFrame(hcat(t_window_total, δ_matrix), vcat("t", delta_names))

    # δ values in relation to the COI
    δ_COI_matrix = zeros(Float64, length(t_window_total), length(δ_COIt))
    aux_count = 0
    for (gen_id, values) in δ_COIt
        aux_count += 1
        δ_COI_matrix[:,aux_count] = rad2deg.(values)
    end
    df_δ_COI = DataFrame(hcat(t_window_total, δ_COI_matrix), vcat("t", gen_names))
    
    # Absolute values of Δω
    Δω_matrix = zeros(Float64, length(t_window_total), length(Δωt)+1)
    aux_count = 0
    for (gen_id, values) in Δωt
        aux_count += 1
        Δω_matrix[:,aux_count] = values
    end
    Δω_matrix[:,end] = ΔωCOIt
    df_Δω = DataFrame(hcat(t_window_total, Δω_matrix), vcat("t", omega_names))

    # Δω values in relation to the COI
    Δω_COI_matrix = zeros(Float64, length(t_window_total), length(Δω_COIt))
    aux_count = 0
    for (gen_id, values) in Δω_COIt
        aux_count += 1
        Δω_COI_matrix[:,aux_count] = values
    end
    df_Δω_COI = DataFrame(hcat(t_window_total, Δω_COI_matrix), vcat("t", gen_names))

    # Frequency values
    f_matrix = zeros(Float64, length(t_window_total), length(Δωt)+1)
    aux_count = 0
    for (gen_id, values) in Δωt
        aux_count += 1
        f_matrix[:,aux_count] = f_syn .* (1.0 .+ values)
    end
    f_matrix[:,end] = f_syn .* (1.0 .+ ΔωCOIt)
    df_f = DataFrame(hcat(t_window_total, f_matrix), vcat("t", omega_names))

    # Frequency values in relation to the COI
    f_COI_matrix = zeros(Float64, length(t_window_total), length(Δω_COIt))
    aux_count = 0
    for (gen_id, values) in Δω_COIt
        aux_count += 1
        f_COI_matrix[:,aux_count] = f_syn .* (1.0 .+ values)
    end
    df_f_COI = DataFrame(hcat(t_window_total, f_COI_matrix), vcat("t", gen_names))
    
    # RoCoF values
    RoCoF_matrix = zeros(Float64, length(time_RoCoF), length(RoCoF)+1)
    aux_count = 0
    for (gen_id, values) in RoCoF
        aux_count += 1
        RoCoF_matrix[:,aux_count] = values
    end
    RoCoF_matrix[:,end] = RoCoFCOI
    df_RoCoF = DataFrame(hcat(time_RoCoF, RoCoF_matrix), vcat("t", omega_names))

    # Vke values
    Vke_matrix = zeros(Float64, length(t_window_total), length(Vke))
    aux_count = 0
    for (gen_id, values) in Vke
        aux_count += 1
        Vke_matrix[:,aux_count] = values
    end
    df_Vke = DataFrame(hcat(t_window_total, Vke_matrix), vcat("t", gen_names))

    # Vpe values
    Vpe_matrix = zeros(Float64, length(t_window_total), length(Vpe))
    aux_count = 0
    for (gen_id, values) in Vpe
        aux_count += 1
        Vpe_matrix[:,aux_count] = values
    end
    df_Vpe = DataFrame(hcat(t_window_total, Vpe_matrix), vcat("t", gen_names))

    # ========================================================
    # Saving the CSVs
    # ========================================================
    cd(joinpath(path_folder_results, "Transient_Stability\\CSV"))
    CSV.write("electrical_power.csv", df_Pe; delim=';')
    CSV.write("accelerating_power.csv", df_Pacc; delim=';')
    CSV.write("accelerating_power_COI.csv", df_Pacc_COI; delim=';')
    CSV.write("angle_abs.csv", df_δ; delim=';')
    CSV.write("angle_rel_COI.csv", df_δ_COI; delim=';')
    CSV.write("speed_dev.csv", df_Δω; delim=';')
    CSV.write("speed_dev_rel_COI.csv", df_Δω_COI; delim=';')
    CSV.write("frequency.csv", df_f; delim=';')
    CSV.write("frequency_rel_COI.csv", df_f_COI; delim=';')
    CSV.write("RoCoF.csv", df_RoCoF; delim=';')
    CSV.write("kinetic_energy.csv", df_Vke; delim=';')
    CSV.write("potential_energy.csv", df_Vpe; delim=';')

    cd(current_path_folder)

    println("Results of the dynamic model successfully saved as CSV files in: ", joinpath(path_folder_results, "Transient_Stability\\CSV"))

end

