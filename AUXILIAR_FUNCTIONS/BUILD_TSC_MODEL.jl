# Functions created to Build the Transient Stability Constraints into the Model

# ===================================================================================
#                          DEFINE INITIAL VARIABLES PRE-FAULT
# ===================================================================================
# Function used to define E, δ and Pm variables of each active generator
function Define_Dyn_Var_EδPm!(model::Model,
    DGEN::DataFrame,
    DGEN_DYN::DataFrame,
    nGEN::Int64,
    base_MVA::Float64,
    V::OrderedDict{Int, JuMP.VariableRef},
    θ::OrderedDict{Int, JuMP.VariableRef},
    P_g::OrderedDict{Int, JuMP.VariableRef},
    Q_g::OrderedDict{Int, JuMP.VariableRef}
    )

    # ======================================================================================
    # Setting the variables of initial internal voltage magnitude and angle of the generator
    # ======================================================================================

    E = OrderedDict{Int, JuMP.VariableRef}()  # Initialize the dictionary to save the variables of internal voltage magnitude
    δ = OrderedDict{Int, JuMP.VariableRef}()  # Initialize the dictionary to save the variables of internal voltage angle
    Pm = OrderedDict{Int, JuMP.VariableRef}() # Initialize the dictionary to save the variables of internal voltage angle

    # Loop to create the variables
    for i in 1:nGEN
        if DGEN.g_status[i] == 1
            E[i] = @variable(model,  # Define the variable related to the internal voltage magnitude
                lower_bound = 0.5,   # Set Lower Bounds
                upper_bound = 2.0,   # Set Upper Bounds
                base_name = "E[$i]", # Set a name
                start = 1.0          # Set a flat start
            )
            δ[i] = @variable(model,  # Define the variable related to the initial internal angle
                lower_bound = -π,    # Set Lower Bounds
                upper_bound = π,     # Set Upper Bounds
                base_name = "δ[$i]", # Set a name
                start = 0.0          # Set a flat start
            )

            Pm[i] = @variable(model,                     # Define the variable related to the initial mechanical power
                lower_bound = DGEN.pg_min[i] / base_MVA, # Set Lower Bounds
                upper_bound = DGEN.pg_max[i] / base_MVA, # Set Upper Bounds
                base_name = "Pm[$i]",                    # Set a name
                start = DGEN.pg_min[i] / base_MVA        # Set a flat start
            )
        end
    end

    # ======================================================================================
    #           Setting the constraints of initial active and reactive
    #               power to define the initial value of E and δ
    # ======================================================================================

    eq_const_P_init = OrderedDict{Int, JuMP.ConstraintRef}() # Initialize the dictionary to save the constraints of active power generated in the pre-fault stage
    eq_const_Q_init = OrderedDict{Int, JuMP.ConstraintRef}() # Initialize the dictionary to save the constraints of reactive power generated in the pre-fault stage
    eq_const_Pm     = OrderedDict{Int, JuMP.ConstraintRef}() # Initialize the dictionary to save the equality constraint of mechanical power in the pre-fault stage

    # Loop to create the equality constraints
    for i in 1:nGEN
        if DGEN.g_status[i] == 1
            bus = DGEN.bus[i]         # ID of the bus in which the generator is connected
            Xd_tr = DGEN_DYN.Xd_tr[i] # Transient reactance of the generator

            # Constraint for active power
            eq_const_P_init[i] = JuMP.@constraint(model,
                (E[i] * V[bus] * sin(δ[i] - θ[bus])) / Xd_tr == P_g[i]
            )

            # Constraint for reactive power
            eq_const_Q_init[i] = JuMP.@constraint(model,
                (E[i] * V[bus] * cos(δ[i] - θ[bus])) / Xd_tr
                - (V[bus]^2) / Xd_tr == Q_g[i]
            )

            # Constraint for mechanical power
            eq_const_Pm[i] = JuMP.@constraint(model,
            Pm[i] == P_g[i]
            )
        end
    end

    return model, E, δ, Pm, eq_const_P_init, eq_const_Q_init, eq_const_Pm
end

# ===================================================================================
#                                  FAULT
# ===================================================================================

# Function used to define the variables that vary in time in the fault period
function Def_Dyn_Fault_All!(
    model::Model,
    DGEN::DataFrame,
    DGEN_DYN::DataFrame,
    nGEN::Int64,
    t_window_fault::Vector{Float64},
    δ_tol::Tuple{Float64, Float64},
    Yred_fault::Matrix,
    E::OrderedDict{Int64, VariableRef},
    Pm::OrderedDict{Int64, VariableRef},
    δ::OrderedDict{Int64, VariableRef},
    Δω_0::Float64,
    P_g::OrderedDict{Int64, VariableRef},
    ω_syn::Float64,
    Δt::Float64
	)

    # ======================================================================================
    #               Setting the variables used in the fault period
    # ======================================================================================

    Pe_tf = OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}() # Initialize the dictionary to save the variables electrical power for each period of time (fault period)
    δ_tf  = OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}() # Initialize the dictionary to save the variables internal angle for each period of time (fault period)
    Δω_tf = OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}() # Initialize the dictionary to save the variables speed deviation for each period of time (fault period)

    # Terms that will be summed to calculate the dynamics of the COI
    terms_δCOI_tf_dict  = OrderedDict{Int, OrderedDict{Int, JuMP.GenericAffExpr{Float64, JuMP.VariableRef}}}() # Dictionary
    terms_ΔωCOI_tf_dict = OrderedDict{Int, OrderedDict{Int, JuMP.GenericAffExpr{Float64, JuMP.VariableRef}}}() # Dictionary

    eq_const_Pe_tf_init = OrderedDict{Int, JuMP.ConstraintRef}()
    eq_const_δ_tf_init  = OrderedDict{Int, JuMP.ConstraintRef}()
    eq_const_Δω_tf_init = OrderedDict{Int, JuMP.ConstraintRef}()

    # Loop to create the variables
    for i in 1:nGEN
        if DGEN.g_status[i] == 1
            # Initialize inner dict for this generator
            Pe_tf[i]     = OrderedDict{Int, JuMP.VariableRef}()
            δ_tf[i]      = OrderedDict{Int, JuMP.VariableRef}()
            Δω_tf[i]     = OrderedDict{Int, JuMP.VariableRef}()

            terms_δCOI_tf_dict[i]  = OrderedDict{Int, JuMP.GenericAffExpr{Float64, JuMP.VariableRef}}()
            terms_ΔωCOI_tf_dict[i] = OrderedDict{Int, JuMP.GenericAffExpr{Float64, JuMP.VariableRef}}()

            for t in eachindex(t_window_fault)
                Pe_tf[i][t] = JuMP.@variable(model, # Define the variable related to the electrical power during fault
                base_name = "Pe_tf[$i, $t]"          # Set the name
                )

                δ_tf[i][t] = JuMP.@variable(model, # Define the variable related to the internal angle during fault
                base_name = "δ_tf[$i, $t]"         # Set the name
                )

                Δω_tf[i][t] = JuMP.@variable(model, # Define the variable related to the speed deviation during fault
                base_name = "Δω_tf[$i, $t]"         # Set the name
                )

                terms_δCOI_tf_dict[i][t]  = DGEN_DYN.H[i] * δ_tf[i][t]
                terms_ΔωCOI_tf_dict[i][t] = DGEN_DYN.H[i] * Δω_tf[i][t]

            end
        end
    end

    # Gives the expressions to calculate the angle and speed deviation of the COI
    expr_δCOI_per_time  = OrderedDict(t => sum(inner_dict[t] for inner_dict in values(terms_δCOI_tf_dict))  for t in eachindex(t_window_fault))
    expr_ΔωCOI_per_time = OrderedDict(t => sum(inner_dict[t] for inner_dict in values(terms_ΔωCOI_tf_dict)) for t in eachindex(t_window_fault))

    # Define the variables and equality constraints of the COI
    model, δCOI_tf, ΔωCOI_tf, eq_const_δCOI_tf, eq_const_ΔωCOI_tf = Def_Dyn_Fault_COI!(model, DGEN_DYN, nGEN, t_window_fault, expr_δCOI_per_time, expr_ΔωCOI_per_time)

    # Define the inequality constraints of angle and speed deviation related to their limits
    model, ineq_const_δ_COI_tf = Def_Dyn_Fault_rel_COI!(model, δ_tf, Δω_tf, δCOI_tf, ΔωCOI_tf, nGEN, t_window_fault, δ_tol)

    # Define the equality constraints for Pe in the fault period
    model, eq_const_Pe_tf = Def_Dyn_Fault_eqconst_Pe!(model, DGEN, nGEN, E, Pe_tf, δ_tf, t_window_fault, Yred_fault)

    # Define the equality constraints of the swing equation for δ and Δω during the fault period
    model, eq_const_δ_tf, eq_const_Δω_tf = Def_Dyn_Fault_Swing!(model, DGEN, DGEN_DYN, nGEN, Pm, Pe_tf, δ_tf, Δω_tf, δ, Δω_0, P_g, t_window_fault, ω_syn, Δt)

    return model, Pe_tf, δ_tf, Δω_tf, δCOI_tf, ΔωCOI_tf, eq_const_Pe_tf_init, eq_const_δ_tf_init, eq_const_Δω_tf_init, eq_const_δCOI_tf, eq_const_ΔωCOI_tf, ineq_const_δ_COI_tf, eq_const_Pe_tf, eq_const_δ_tf, eq_const_Δω_tf
end

# Function that define the variables of the COI as well as its equality constraints across the fault period
function Def_Dyn_Fault_COI!(model::Model,
    DGEN_DYN::DataFrame,
    nGEN::Int64,
    t_window_fault::Vector{Float64},
    aux_expr_δCOI_per_time::OrderedDict{Int64, AffExpr},
    aux_expr_ΔωCOI_per_time::OrderedDict{Int64, AffExpr}
    )

    H_COI = sum(DGEN_DYN.H[i] for i in 1:nGEN if DGEN.g_status[i] == 1) # Inertia of the COI

    δCOI_tf  = OrderedDict{Int, JuMP.VariableRef}() # Initialize the dictionary to save the variables angle of the COI for each period of time (fault period)
    ΔωCOI_tf = OrderedDict{Int, JuMP.VariableRef}() # Initialize the dictionary to save the variables speed deviation of the COI for each period of time (fault period)

    eq_const_δCOI_tf  = OrderedDict{Int, JuMP.ConstraintRef}()
    eq_const_ΔωCOI_tf = OrderedDict{Int, JuMP.ConstraintRef}()

    for t in eachindex(t_window_fault)
        δCOI_tf[t]  = JuMP.@variable(model, base_name = "δCOI_tf[$t]")
        ΔωCOI_tf[t] = JuMP.@variable(model, base_name = "ΔωCOI_tf[$t]")
        
        eq_const_δCOI_tf[t]  = JuMP.@constraint(model, δCOI_tf[t]  == aux_expr_δCOI_per_time[t]  / H_COI)
        eq_const_ΔωCOI_tf[t] = JuMP.@constraint(model, ΔωCOI_tf[t] == aux_expr_ΔωCOI_per_time[t] / H_COI)

    end

    return model, δCOI_tf, ΔωCOI_tf, eq_const_δCOI_tf, eq_const_ΔωCOI_tf

end

# Function that define the constraints of the angle of generators in relation to the COI variables across the fault period
function Def_Dyn_Fault_rel_COI!(model::Model,
    δ_tf::OrderedDict{Int64, OrderedDict{Int64, VariableRef}},
    Δω_tf::OrderedDict{Int64, OrderedDict{Int64, VariableRef}},
    δCOI_tf::OrderedDict{Int, JuMP.VariableRef},
    ΔωCOI_tf::OrderedDict{Int, JuMP.VariableRef},
    nGEN::Int64,
    t_window_fault::Vector{Float64},
    δ_tol::Tuple{Float64, Float64}
    )

    ineq_const_δ_COI_tf  = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}()

    # Loop to create the variables
    for i in 1:nGEN
        if DGEN.g_status[i] == 1
            # Initialize inner dict for this generator
            ineq_const_δ_COI_tf[i]  = OrderedDict{Int, JuMP.ConstraintRef}()

            for t in eachindex(t_window_fault)
                ineq_const_δ_COI_tf[i][t] = JuMP.@constraint(model,
                δ_tol[1] <= δ_tf[i][t] - δCOI_tf[t] <= δ_tol[2]
                )
            end
        end
    end

    return model, ineq_const_δ_COI_tf


end

# Function that create the constraints of the elecrical power across the fault period
function Def_Dyn_Fault_eqconst_Pe!(model::Model,
    DGEN::DataFrame,
    nGEN::Int64,
    E::OrderedDict{Int64, VariableRef},
    Pe_tf::OrderedDict{Int64, OrderedDict{Int64, VariableRef}},
    δ_tf::OrderedDict{Int64, OrderedDict{Int64, VariableRef}},
    t_window_fault::Vector{Float64},
    Yred_fault::Matrix
    )
    
    eq_const_Pe_tf = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}() # Dictionary of equality constraints

    terms_Pe_tf_dict = OrderedDict{Int, OrderedDict{Int, JuMP.NonlinearExpr}}() # Dictionary of Terms that will be summed to calculate Pe

    active_gen = findall(isone, DGEN.g_status) # Vector with the indices of the active generators

    for (i, id_gen1) in enumerate(active_gen)

        eq_const_Pe_tf[id_gen1] = OrderedDict{Int, JuMP.ConstraintRef}() # Dictionary of equality constraints

        terms_Pe_tf_dict[id_gen1] = OrderedDict{Int, JuMP.NonlinearExpr}() # Dictionary of terms to form a nonlinear expression

        for t in eachindex(t_window_fault) # Loop over the time
            terms_to_sum = JuMP.NonlinearExpr[] # List of expressions

            for (j, id_gen2) in enumerate(active_gen) #in eachindex(active_gen) # Loop over active generators
                G = real(Yred_fault[i,j]) # Conductance of the reduced admittance matrix
                B = imag(Yred_fault[i,j]) # Susceptance of the reduced admittance matrix

                expression_Pe = E[id_gen2] * (G*cos(δ_tf[id_gen1][t] - δ_tf[id_gen2][t]) + B*sin(δ_tf[id_gen1][t] - δ_tf[id_gen2][t]))
                push!(terms_to_sum, expression_Pe)
            end
            terms_Pe_tf_dict[id_gen1][t] = E[id_gen1] * sum(terms_to_sum)
            eq_const_Pe_tf[id_gen1][t] = JuMP.@constraint(model, Pe_tf[id_gen1][t] == terms_Pe_tf_dict[id_gen1][t]) # Equality constraint for the electrical power
        end
    end

    return model, eq_const_Pe_tf

end

# Function that model the swing equation using Trapezoidal Approximation
function Def_Dyn_Fault_Swing!(model::Model,
    DGEN::DataFrame,
    DGEN_DYN::DataFrame,
    nGEN::Int64,
    Pm::OrderedDict{Int64, VariableRef},
    Pe_tf::OrderedDict{Int64, OrderedDict{Int64, VariableRef}},
    δ_tf::OrderedDict{Int64, OrderedDict{Int64, VariableRef}},
    Δω_tf::OrderedDict{Int64, OrderedDict{Int64, VariableRef}},
    δ::OrderedDict{Int64, VariableRef},
    Δω_0::Float64,
    P_g::OrderedDict{Int64, VariableRef},
    t_window_fault::Vector{Float64},
    ω_syn::Float64,
    Δt::Float64
    )

    eq_const_δ_tf  = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}() # Dictionary of equality constraints
    eq_const_Δω_tf = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}() # Dictionary of equality constraints

    for i in 1:nGEN
        if DGEN.g_status[i] == 1

            eq_const_δ_tf[i] = OrderedDict{Int, JuMP.ConstraintRef}()
            eq_const_Δω_tf[i] = OrderedDict{Int, JuMP.ConstraintRef}()

            H = DGEN_DYN.H[i]
            D = DGEN_DYN.D[i]

            for t in eachindex(t_window_fault)
                # ====================
                # Trapezoidal Method
                # ====================
                if t == 1
                    eq_const_δ_tf[i][t] = JuMP.@constraint(model, δ_tf[i][t] - δ[i] - (ω_syn * (Δt/2) * (Δω_tf[i][t] + Δω_0)) == 0.0)

                    eq_const_Δω_tf[i][t] = JuMP.@constraint(model, ((1.0 + ((D*Δt) / (4*H))) * Δω_tf[i][t]) - ((1.0 - ((D*Δt) / (4*H))) * Δω_0) - (Δt / (4*H))*(2*Pm[i] - Pe_tf[i][t] - P_g[i])  == 0.0)
                else
                    eq_const_δ_tf[i][t] = JuMP.@constraint(model, δ_tf[i][t] - δ_tf[i][t-1] - (ω_syn * (Δt/2) * (Δω_tf[i][t] + Δω_tf[i][t-1])) == 0.0)

                    eq_const_Δω_tf[i][t] = JuMP.@constraint(model, ((1.0 + ((D*Δt) / (4*H))) * Δω_tf[i][t]) - ((1.0 - ((D*Δt) / (4*H))) * Δω_tf[i][t-1]) - (Δt / (4*H))*(2*Pm[i] - Pe_tf[i][t] - Pe_tf[i][t-1])  == 0.0)
                end
            end
        end
    end
    return model, eq_const_δ_tf, eq_const_Δω_tf

end

# ===================================================================================
#                               POST-FAULT
# ===================================================================================

# Function used to define the variables that vary in time in the post-fault period
function Def_Dyn_PostF_All!(model::Model,
    DGEN::DataFrame,
    DGEN_DYN::DataFrame,
    nGEN::Int64,
    t_window_postf::Vector{Float64},
    δ_tol::Tuple{Float64, Float64},
    Yred_postf::Matrix{ComplexF64},
    E::OrderedDict{Int64, VariableRef},
    Pm::OrderedDict{Int64, VariableRef},
    δ_tf::OrderedDict{Int64, OrderedDict{Int64, VariableRef}},
    Δω_tf::OrderedDict{Int64, OrderedDict{Int64, VariableRef}},
    Pe_tf::OrderedDict{Int64, OrderedDict{Int64, VariableRef}},
    ω_syn::Float64,
    Δt::Float64
    )

    # ======================================================================================
    #               Setting the variables used in the post-fault period
    # ======================================================================================

    Pe_tpf = OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}() # Initialize the dictionary to save the variables electrical power for each period of time (post-fault period)
    δ_tpf  = OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}() # Initialize the dictionary to save the variables internal angle for each period of time (post-fault period)
    Δω_tpf = OrderedDict{Int, OrderedDict{Int, JuMP.VariableRef}}() # Initialize the dictionary to save the variables speed deviation for each period of time (post-fault period)

    # Terms that will be summed to calculate the dynamics of the COI
    terms_δCOI_tpf_dict  = OrderedDict{Int, OrderedDict{Int, JuMP.GenericAffExpr{Float64, JuMP.VariableRef}}}() # Dictionary
    terms_ΔωCOI_tpf_dict = OrderedDict{Int, OrderedDict{Int, JuMP.GenericAffExpr{Float64, JuMP.VariableRef}}}() # Dictionary

    eq_const_Pe_tpf_init = OrderedDict{Int, JuMP.ConstraintRef}()
    eq_const_δ_tpf_init  = OrderedDict{Int, JuMP.ConstraintRef}()
    eq_const_Δω_tpf_init = OrderedDict{Int, JuMP.ConstraintRef}()

    # Loop to create the variables
    for i in 1:nGEN
        if DGEN.g_status[i] == 1
            # Initialize inner dict for this generator
            Pe_tpf[i]     = OrderedDict{Int, JuMP.VariableRef}()
            δ_tpf[i]      = OrderedDict{Int, JuMP.VariableRef}()
            Δω_tpf[i]     = OrderedDict{Int, JuMP.VariableRef}()

            terms_δCOI_tpf_dict[i]  = OrderedDict{Int, JuMP.GenericAffExpr{Float64, JuMP.VariableRef}}()
            terms_ΔωCOI_tpf_dict[i] = OrderedDict{Int, JuMP.GenericAffExpr{Float64, JuMP.VariableRef}}()

            for t in eachindex(t_window_postf)
                Pe_tpf[i][t] = JuMP.@variable(model, # Define the variable related to the electrical power during post-fault
                base_name = "Pe_tpf[$i, $t]"          # Set the name
                )

                δ_tpf[i][t] = JuMP.@variable(model, # Define the variable related to the internal angle during post-fault
                base_name = "δ_tpf[$i, $t]"         # Set the name
                )

                Δω_tpf[i][t] = JuMP.@variable(model, # Define the variable related to the speed deviation during post-fault
                base_name = "Δω_tpf[$i, $t]"         # Set the name
                )

                terms_δCOI_tpf_dict[i][t]  = DGEN_DYN.H[i] * δ_tpf[i][t]
                terms_ΔωCOI_tpf_dict[i][t] = DGEN_DYN.H[i] * Δω_tpf[i][t]

            end
        end
    end

    # Gives the expressions to calculate the angle and speed deviation of the COI
    expr_δCOI_per_time  = OrderedDict(t => sum(inner_dict[t] for inner_dict in values(terms_δCOI_tpf_dict))  for t in eachindex(t_window_postf))
    expr_ΔωCOI_per_time = OrderedDict(t => sum(inner_dict[t] for inner_dict in values(terms_ΔωCOI_tpf_dict)) for t in eachindex(t_window_postf))

    # Define the variables and equality constraints of the COI
    model, δCOI_tpf, ΔωCOI_tpf, eq_const_δCOI_tpf, eq_const_ΔωCOI_tpf = Def_Dyn_PostF_COI!(model, DGEN_DYN, nGEN, t_window_postf, expr_δCOI_per_time, expr_ΔωCOI_per_time)

    # Define the inequality constraints of angle and speed deviation related to their limits
    model, ineq_const_δ_COI_tpf = Def_Dyn_PostF_rel_COI!(model, δ_tpf, Δω_tpf, δCOI_tpf, ΔωCOI_tpf, nGEN, t_window_postf, δ_tol)

    # Define the equality constraints for Pe in the post-fault period
    model, eq_const_Pe_tpf = Def_Dyn_PostF_eqconst_Pe!(model, DGEN, nGEN, E, Pe_tpf, δ_tpf, t_window_postf, Yred_postf)

    # Define the equality constraints of the swing equation for δ and Δω during the post-fault period
    model, eq_const_δ_tpf, eq_const_Δω_tpf = Def_Dyn_PostF_Swing!(model, DGEN, DGEN_DYN, nGEN, Pm, Pe_tpf, δ_tpf, Δω_tpf, δ_tf, Δω_tf, Pe_tf, t_window_postf, ω_syn, Δt)

    return model, Pe_tpf, δ_tpf, Δω_tpf, δCOI_tpf, ΔωCOI_tpf, eq_const_Pe_tpf_init, eq_const_δ_tpf_init, eq_const_Δω_tpf_init, eq_const_δCOI_tpf, eq_const_ΔωCOI_tpf, ineq_const_δ_COI_tpf, eq_const_Pe_tpf, eq_const_δ_tpf, eq_const_Δω_tpf
end

# Function that define the variables of the COI as well as its equality constraints across the post-fault period
function Def_Dyn_PostF_COI!(model::Model,
    DGEN_DYN::DataFrame,
    nGEN::Int64,
    t_window_postf::Vector{Float64},
    aux_expr_δCOI_per_time::OrderedDict{Int64, AffExpr},
    aux_expr_ΔωCOI_per_time::OrderedDict{Int64, AffExpr}
    )

    H_COI = sum(DGEN_DYN.H[i] for i in 1:nGEN if DGEN.g_status[i] == 1) # Inertia of the COI

    δCOI_tpf  = OrderedDict{Int, JuMP.VariableRef}() # Initialize the dictionary to save the variables angle of the COI for each period of time (fault period)
    ΔωCOI_tpf = OrderedDict{Int, JuMP.VariableRef}() # Initialize the dictionary to save the variables speed deviation of the COI for each period of time (fault period)

    eq_const_δCOI_tpf  = OrderedDict{Int, JuMP.ConstraintRef}()
    eq_const_ΔωCOI_tpf = OrderedDict{Int, JuMP.ConstraintRef}()

    for t in eachindex(t_window_postf)
        δCOI_tpf[t]  = JuMP.@variable(model, base_name = "δCOI_tpf[$t]",  start = 0.0)
        ΔωCOI_tpf[t] = JuMP.@variable(model, base_name = "ΔωCOI_tpf[$t]", start = 0.0)
        
        eq_const_δCOI_tpf[t]  = JuMP.@constraint(model, δCOI_tpf[t]  == aux_expr_δCOI_per_time[t]  / H_COI)
        eq_const_ΔωCOI_tpf[t] = JuMP.@constraint(model, ΔωCOI_tpf[t] == aux_expr_ΔωCOI_per_time[t] / H_COI)

    end

    return model, δCOI_tpf, ΔωCOI_tpf, eq_const_δCOI_tpf, eq_const_ΔωCOI_tpf

end

# Function that define the constraints of the angle of generators in relation to the COI variables across the post-fault period
function Def_Dyn_PostF_rel_COI!(model::Model,
    δ_tpf::OrderedDict{Int64, OrderedDict{Int64, VariableRef}},
    Δω_tpf::OrderedDict{Int64, OrderedDict{Int64, VariableRef}},
    δCOI_tpf::OrderedDict{Int, JuMP.VariableRef},
    ΔωCOI_tpf::OrderedDict{Int, JuMP.VariableRef},
    nGEN::Int64,
    t_window_postf::Vector{Float64},
    δ_tol::Tuple{Float64, Float64},
    )

    ineq_const_δ_COI_tpf  = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}()

    # Loop to create the variables
    for i in 1:nGEN
        if DGEN.g_status[i] == 1
            # Initialize inner dict for this generator
            ineq_const_δ_COI_tpf[i]  = OrderedDict{Int, JuMP.ConstraintRef}()

            for t in eachindex(t_window_postf)
                ineq_const_δ_COI_tpf[i][t] = JuMP.@constraint(model,
                δ_tol[1] <= δ_tpf[i][t] - δCOI_tpf[t] <= δ_tol[2]
                )
            end
        end
    end

    return model, ineq_const_δ_COI_tpf


end

# Function that create the constraints of the elecrical power across the post-fault period
function Def_Dyn_PostF_eqconst_Pe!(model::Model,
    DGEN::DataFrame,
    nGEN::Int64,
    E::OrderedDict{Int64, VariableRef},
    Pe_tpf::OrderedDict{Int64, OrderedDict{Int64, VariableRef}},
    δ_tpf::OrderedDict{Int64, OrderedDict{Int64, VariableRef}},
    t_window_postf::Vector{Float64},
    Yred_postf::Matrix{ComplexF64}
    )
    
    eq_const_Pe_tpf = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}() # Dictionary of equality constraints

    terms_Pe_tpf_dict = OrderedDict{Int, OrderedDict{Int, JuMP.NonlinearExpr}}() # Dictionary of Terms that will be summed to calculate Pe

    active_gen = findall(isone, DGEN.g_status) # Vector with the indices of the active generators

    for (i, id_gen1) in enumerate(active_gen)

        eq_const_Pe_tpf[id_gen1] = OrderedDict{Int, JuMP.ConstraintRef}() # Dictionary of equality constraints

        terms_Pe_tpf_dict[id_gen1] = OrderedDict{Int, JuMP.NonlinearExpr}() # Dictionary of terms to form a nonlinear expression

        for t in eachindex(t_window_postf) # Loop over the time
            terms_to_sum = JuMP.NonlinearExpr[] # List of expressions

            for (j, id_gen2) in enumerate(active_gen) #in eachindex(active_gen) # Loop over active generators
                G = real(Yred_postf[i,j]) # Conductance of the reduced admittance matrix
                B = imag(Yred_postf[i,j]) # Susceptance of the reduced admittance matrix

                expression_Pe = E[id_gen2] * (G*cos(δ_tpf[id_gen1][t] - δ_tpf[id_gen2][t]) + B*sin(δ_tpf[id_gen1][t] - δ_tpf[id_gen2][t]))
                push!(terms_to_sum, expression_Pe)
            end
            terms_Pe_tpf_dict[id_gen1][t] = E[id_gen1] * sum(terms_to_sum)
            eq_const_Pe_tpf[id_gen1][t] = JuMP.@constraint(model, Pe_tpf[id_gen1][t] == terms_Pe_tpf_dict[id_gen1][t]) # Equality constraint for the electrical power
        end
    end

    return model, eq_const_Pe_tpf

end

# Function that model the swing equation using Trapezoidal Approximation (Post-Fault period)
function Def_Dyn_PostF_Swing!(model::Model,
    DGEN::DataFrame,
    DGEN_DYN::DataFrame,
    nGEN::Int64,
    Pm::OrderedDict{Int64, VariableRef},
    Pe_tpf::OrderedDict{Int64, OrderedDict{Int64, VariableRef}},
    δ_tpf::OrderedDict{Int64, OrderedDict{Int64, VariableRef}},
    Δω_tpf::OrderedDict{Int64, OrderedDict{Int64, VariableRef}},
    δ_tf::OrderedDict{Int64, OrderedDict{Int64, VariableRef}},
    Δω_tf::OrderedDict{Int64, OrderedDict{Int64, VariableRef}},
    Pe_tf::OrderedDict{Int64, OrderedDict{Int64, VariableRef}},
    t_window_postf::Vector{Float64},
    ω_syn::Float64,
    Δt::Float64
    )

    eq_const_δ_tpf  = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}() # Dictionary of equality constraints
    eq_const_Δω_tpf = OrderedDict{Int, OrderedDict{Int, JuMP.ConstraintRef}}() # Dictionary of equality constraints

    for i in 1:nGEN
        if DGEN.g_status[i] == 1

            eq_const_δ_tpf[i] = OrderedDict{Int, JuMP.ConstraintRef}()
            eq_const_Δω_tpf[i] = OrderedDict{Int, JuMP.ConstraintRef}()

            H = DGEN_DYN.H[i]
            D = DGEN_DYN.D[i]
            
            for t in eachindex(t_window_postf)
                # ====================
                # Trapezoidal Method
                # ====================
                if t == 1
                    last_var_δ_tf = last(δ_tf[i])[2]   # last VariableRef
                    last_var_Δω_tf = last(Δω_tf[i])[2]   # last VariableRef
                    last_var_Pe_tf = last(Pe_tf[i])[2]   # last VariableRef

                    eq_const_δ_tpf[i][t] = JuMP.@constraint(model, δ_tpf[i][t] - last_var_δ_tf - (ω_syn * (Δt/2) * (Δω_tpf[i][t] + last_var_Δω_tf)) == 0.0)

                    eq_const_Δω_tpf[i][t] = JuMP.@constraint(model, ((1.0 + ((D*Δt) / (4*H))) * Δω_tpf[i][t]) - ((1.0 - ((D*Δt) / (4*H))) * last_var_Δω_tf) - (Δt / (4*H))*(2*Pm[i] - Pe_tpf[i][t] - last_var_Pe_tf)  == 0.0)
                else
                    eq_const_δ_tpf[i][t] = JuMP.@constraint(model, δ_tpf[i][t] - δ_tpf[i][t-1] - (ω_syn * (Δt/2) * (Δω_tpf[i][t] + Δω_tpf[i][t-1])) == 0.0)

                    eq_const_Δω_tpf[i][t] = JuMP.@constraint(model, ((1.0 + ((D*Δt) / (4*H))) * Δω_tpf[i][t]) - ((1.0 - ((D*Δt) / (4*H))) * Δω_tpf[i][t-1]) - (Δt / (4*H))*(2*Pm[i] - Pe_tpf[i][t] - Pe_tpf[i][t-1])  == 0.0)
                end
            end
        end
    end
    return model, eq_const_δ_tpf, eq_const_Δω_tpf

end


