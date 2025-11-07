cd(dirname(@__FILE__));

#=
CODE FOR SOLVING THE AC POWER FLOW WITHOUT RELAXATIONS
AND ADDING TRANSIENT STABILITY CONSTRAINTS

Author:      Alex Junior da Cunha Coelho
Supervisors: Luis Badesa Bernardo and Araceli Hernandez Bayo
Affiliation: Technical University of Madrid
October 2025
=#

#---------------------------
# INCLUDE THE PACKAGES USED
#---------------------------
# Packages related to linear algebra
using LinearAlgebra, SparseArrays

# Packages related to treatement of data
using Dates, NumericIO, DataFrames, Printf, CSV, DataStructures

# Packages related to the optimization
using JuMP, Ipopt, MathOptInterface

# Packages for plotting
using Plots, LaTeXStrings, Measures

# Packages for numerical integration
using Trapz

#---------------------------------
# INCLUDE AUXILIAR FUNCTION FILES
#--------------------------------
include("./AUXILIAR_FUNCTIONS//AF_CLEAN_TERMINAL.jl")         # Auxiliar function to clean the terminal
include("./AUXILIAR_FUNCTIONS//AF_DERIVATIVE_TECHINIQUES.jl") # Auxiliar functions to calculate numerical derivatives
include("./AUXILIAR_FUNCTIONS//AF_READ_DATA.jl")              # Auxiliar functions used to read input data
include("./AUXILIAR_FUNCTIONS//AF_MANAGEMENT.jl")             # Auxiliar functions to calculate AC power flow and manage some data
include("./AUXILIAR_FUNCTIONS//AF_YBUS.jl")                   # Auxiliar function to create Ybus
include("./AUXILIAR_FUNCTIONS//BUILD_ACOPF_MODEL.jl")         # Auxiliar function to create the AC OPF model for optimization
include("./AUXILIAR_FUNCTIONS//BUILD_TSC_MODEL.jl")           # Auxiliar function to create to include the Transient Stability Constraints into the optimization model
include("./AUXILIAR_FUNCTIONS//AF_SAVE_OUTPUT.jl")            # Auxiliar function to save the output results

Clean_Terminal() # Clean the terminal

#-----------------------------------------
# Generate a folder to export the results
#-----------------------------------------
current_path_folder = pwd()                                            # Directory of the current folder
name_path_results   = "RESULTS"                                        # Name of the folder to save the results (it must be created in advance)
path_folder_results = joinpath(current_path_folder, name_path_results) # Results directory
cd(current_path_folder)                                                # Load the current folder

#=
----------------------------------------------
 Relevant input variables to solve the TSC-OPF
----------------------------------------------
** Select a system from the options below: **
9bus
=#

case     = "9bus" # Case under study (folder name)
base_MVA = 100.0  # Base Power [MVA]
δ_tol  = 100 # Module of maximum tolerance for internal angle [degrees]
δ_tol  = (deg2rad(-δ_tol), deg2rad(δ_tol)) # Tupple with the min and max tolerance for internal angle

# System Parameters
f_syn = 50.0        # Synchronous frequency in [Hz]
ω_syn = 2π * f_syn  # Synchronous speed in     [rad/s]
Δω_0  = 0.0         # Initial speed deviation  [p.u.]

# # ----------- Contingency Number 1 -------------- #
# bus_fault = 4   # Faulted bus (3Φ short circuit)
# circ_trip = [4] # Circuit(s) that trip(s) and clear(s) the fault (line 4-5)

# ----------- Contingency Number 2 -------------- #
bus_fault = 7   # Faulted bus (3Φ short circuit)
circ_trip = [6] # Circuit(s) that trip(s) and clear(s) the fault (line 7-5)


# Variables related to the simulation time
t_start_sim    = 0.0    # Time when the simulation is started           [seconds] 
t_end_sim      = 5.0    # Time to stop the simulation                   [seconds]
t_step         = 0.01   # Time step to solve the differential equations [seconds]
t_start_fault  = 0.01   # Time when the fault is started                [seconds]
clearing_time  = 0.30   # Time the fault lasts until clearing           [seconds] 
t_clear_fault  = clearing_time + t_start_fault   # Time to clear the fault [seconds] 
t_window_fault = collect(t_start_fault:t_step:t_clear_fault) # Time window - During-fault period
t_window_postf = collect(t_clear_fault:t_step:t_end_sim)     # Time window - Post-fault period
t_window_total = vcat(t_window_fault, t_window_postf)        # Time window - Total simulation period

if t_end_sim <= t_clear_fault
    throw(ArgumentError("t_end_sim is lower than t_clear_fault; PLEASE CORRECT."))
end

println("--------------------------------------------------------------------------------------------------------------------------------------")
Print_Input_Parameters(case, base_MVA, δ_tol, f_syn, bus_fault, circ_trip, t_start_sim, t_end_sim, t_step, t_start_fault, clearing_time, t_clear_fault, current_path_folder, path_folder_results)

#------------------------------------------
# Call the function to read the input data
#------------------------------------------
input_data_path_folder = joinpath(current_path_folder, "INPUT_DATA", case) # Folder name where the input data is located
# Get the structs with data related to buses, generators and circuits
DBUS, DGEN, DGEN_DYN, DCIR, bus_mapping, reverse_bus_mapping = Read_Input_Data(input_data_path_folder)

# Variable to multiply the power demanded by the loads
load_factor = 1.5
DBUS.p_d = load_factor .* (DBUS.p_d)
DBUS.q_d = load_factor .* (DBUS.q_d)
# ---------------------------------------------------

nBUS = length(DBUS.bus)      # Number of buses in the system
nGEN = length(DGEN.id)       # Number of generators in the system
if nGEN != length(DGEN_DYN.id) throw(ArgumentError("Dynamic and static data of generators do not match in length.")) end
nCIR = length(DCIR.from_bus) # Number of circuits in the system
cd(current_path_folder)      # Load the current folder

# ------------------
# Some sanity checks
# ------------------
if !(haskey(bus_mapping, bus_fault))
    throw(ArgumentError("The ID of the faulted bus does not exist in the CSV of input data."))
end

if maximum(circ_trip) > nCIR
    throw(ArgumentError("The IDs of the faulted circuits are greater than the ids in DCIR."))
end

#-------------------------------------------------------------------------------------------------
# Associates the buses with the generators and circuits connected to it, as well as adjacent buses
#-------------------------------------------------------------------------------------------------
bus_gen_circ_dict, bus_gen_circ_dict_ON = Manage_Bus_Gen_Circ(DBUS, DGEN, DCIR) 

# -----------------------------------
# Calculate the admittance matrices
# -----------------------------------

# Calculate the admittance matrix for the pre-fault period
Ybus_pref = Calculate_Ybus(DBUS, DCIR, nBUS, nCIR, base_MVA) # Admittance matrix for the pre-fault stage

# Calculate the new admittance matrix for the fault period and its reduced version 
Ybus_fault = Calculate_Ybus_fault(deepcopy(Ybus_pref), DBUS, DGEN, DGEN_DYN, nBUS, nGEN, base_MVA, bus_fault)                               # Admittance matrix augmented for the fault
Yred_fault = (Ybus_fault[nBUS+1:end, nBUS+1:end]) - (Ybus_fault[nBUS+1:end, 1:nBUS]) * inv(Matrix(Ybus_fault[1:nBUS, 1:nBUS])) * Ybus_fault[1:nBUS, nBUS+1:end] # Reduced admittance matrix for the fault

# Calculate the new admittance matrix for the post-fault period and its reduced version 
l_status_postf = deepcopy(DCIR.l_status)
l_status_postf[circ_trip] .= 0
Ybus_postf = Calculate_Ybus_postf(DBUS, DCIR, DGEN, DGEN_DYN, nBUS, nCIR, nGEN, base_MVA, l_status_postf)                                     # Admittance matrix augmented for the post-fault
Yred_postf = (Ybus_postf[nBUS+1:end, nBUS+1:end]) - ((Ybus_postf[nBUS+1:end, 1:nBUS]) * inv(Matrix(Ybus_postf[1:nBUS, 1:nBUS])) * Ybus_postf[1:nBUS, nBUS+1:end]) # Reduced admittance matrix for the post-fault


cd(joinpath(path_folder_results,"Admittance_Matrices"))
df_Ybus_pref = DataFrame(Matrix(Ybus_pref), :auto)        # Convert the admittance matrix of the pre-fault into a DataFrame to save it
CSV.write("df_Ybus_pref.csv", df_Ybus_pref; delim=';')    # Save the admittance matrix of the pre-fault in a CSV file

df_Ybus_fault = DataFrame(Matrix(Ybus_fault), :auto)      # Convert the admittance matrix of the fault into a DataFrame to save it
CSV.write("df_Ybus_fault.csv", df_Ybus_fault; delim=';')  # Save the admittance matrix of the fault in a CSV file

df_Yred_fault = DataFrame(Matrix(Yred_fault), :auto)      # Convert the reduced admittance matrix of the fault into a DataFrame to save it
CSV.write("df_Yred_fault.csv", df_Yred_fault; delim=';')  # Save the reduced admittance matrix of the fault in a CSV file

df_Ybus_postf = DataFrame(Matrix(Ybus_postf), :auto)      # Convert the admittance matrix of the post-fault into a DataFrame to save it
CSV.write("df_Ybus_postf.csv", df_Ybus_postf; delim=';')  # Save the admittance matrix of the post-fault in a CSV file

df_Yred_postf = DataFrame(Matrix(Yred_postf), :auto)      # Convert the reduced admittance matrix of the post-fault into a DataFrame to save it
CSV.write("df_Yred_postf.csv", df_Yred_postf; delim=';')  # Save the reduced admittance matrix of the post-fault in a CSV file

println("Admittance matrices successfully saved in: ", joinpath(path_folder_results,"Admittance_Matrices"))
println("--------------------------------------------------------------------------------------------------------------------------------------")
cd(current_path_folder)

# ------------------
# Some sanity checks
# ------------------
active_gen = findall(isone, DGEN.g_status)
if length(active_gen) != size(Yred_fault)[1] && length(active_gen) != size(Yred_fault)[2]
    throw(ArgumentError("The size of the fault reduced matrix must be equal the number of active generators."))
end
if length(active_gen) != size(Yred_postf)[1] && length(active_gen) != size(Yred_postf)[2]
    throw(ArgumentError("The size of the post-fault reduced matrix must be equal the number of active generators."))
end
# ########################################################################################
#                                 STARTS OPTIMIZATION PROCESS 

#-----------------------------------
# Optimization model -> Setup
#-----------------------------------
optimizer = Ipopt.Optimizer
model = JuMP.Model(optimizer)
JuMP.set_optimizer_attribute(model, "tol", 1e-9)
JuMP.set_optimizer_attribute(model, "acceptable_tol", 1e-9)
JuMP.set_optimizer_attribute(model, "print_level", 5)               # Verbosity (0–12, default = 5)
JuMP.set_optimizer_attribute(model, "output_file", "ipopt_log.txt") # Prints to console instead of a file
JuMP.set_optimizer_attribute(model, "max_iter", 5000)
JuMP.set_silent(model)

#------------------------------
# Build the Optimization Model
#------------------------------
time_to_build_model = time() # Start the timer to build the Optimization Model

model, V, θ, P_g, Q_g, P_ik, Q_ik, P_ki, Q_ki, eq_const_angle_sw, eq_const_p_balance, eq_const_q_balance, 
eq_const_p_ik, eq_const_q_ik, eq_const_p_ki, eq_const_q_ki, ineq_const_s_ik, ineq_const_s_ki, 
ineq_const_diff_ang = Make_ACOPF_Model!(model, DBUS, DGEN, DCIR, bus_gen_circ_dict_ON, base_MVA, nBUS, nGEN, nCIR)

# =====================================================================================
# Adding transient stability constraints

# Add variables of the pre-fault condition
model, E, δ, Pm, eq_const_P_init, eq_const_Q_init, eq_const_Pm = Define_Dyn_Var_EδPm!(model, DGEN, DGEN_DYN, nGEN, base_MVA, V, θ, P_g, Q_g)

# ----------------------------------------------------------------------------------------------------
# Add variables and constraints of the fault period
model, Pe_tf, δ_tf, Δω_tf, δCOI_tf, ΔωCOI_tf, eq_const_Pe_tf_init, eq_const_δ_tf_init, eq_const_Δω_tf_init, eq_const_δCOI_tf, eq_const_ΔωCOI_tf, 
ineq_const_δ_COI_tf, eq_const_Pe_tf, eq_const_δ_tf, eq_const_Δω_tf = Def_Dyn_Fault_All!(model, DGEN, DGEN_DYN, nGEN, 
t_window_fault, δ_tol, Yred_fault, E, Pm, δ, Δω_0, P_g, ω_syn, t_step)

# ----------------------------------------------------------------------------------------------------
# Add variables and constraints of the post-fault period
model, Pe_tpf, δ_tpf, Δω_tpf, δCOI_tpf, ΔωCOI_tpf, eq_const_Pe_tpf_init, eq_const_δ_tpf_init, eq_const_Δω_tpf_init, eq_const_δCOI_tpf, 
eq_const_ΔωCOI_tpf, ineq_const_δ_COI_tpf, eq_const_Pe_tpf, eq_const_δ_tpf, eq_const_Δω_tpf = Def_Dyn_PostF_All!(model, 
DGEN, DGEN_DYN, nGEN, t_window_postf, δ_tol, Yred_postf, E, Pm, δ_tf, Δω_tf, Pe_tf, ω_syn, t_step)

time_to_build_model = time() - time_to_build_model # End the timer to build the Optimization Model
println("\nTime to build the model: $time_to_build_model sec\n")

# =====================================================================================

#-------------------------------------------------------------------------------------
#                         SAVE MODEL SUMMARY AND DETAILS
#-------------------------------------------------------------------------------------
println("--------------------------------------------------------------------------------------------------------------------------------------")
Export_ACOPF_Model(model, V, θ, P_g, Q_g, P_ik, Q_ik, P_ki, Q_ki, eq_const_angle_sw, eq_const_p_balance, eq_const_q_balance, 
eq_const_p_ik, eq_const_q_ik, eq_const_p_ki, eq_const_q_ki, ineq_const_s_ik, ineq_const_s_ki, 
ineq_const_diff_ang, current_path_folder, path_folder_results)

Export_Dyn_Model(model, E, δ, Pm, eq_const_Pm, Pe_tf, δ_tf, Δω_tf, δCOI_tf, ΔωCOI_tf, eq_const_δCOI_tf, eq_const_ΔωCOI_tf, eq_const_Pe_tf, eq_const_δ_tf,
eq_const_Δω_tf, ineq_const_δ_COI_tf, Pe_tpf, δ_tpf, Δω_tpf, δCOI_tpf, ΔωCOI_tpf, eq_const_δCOI_tpf, eq_const_ΔωCOI_tpf, eq_const_Pe_tpf,
eq_const_δ_tpf, eq_const_Δω_tpf, ineq_const_δ_COI_tpf, current_path_folder, path_folder_results)
println("--------------------------------------------------------------------------------------------------------------------------------------")

# ---------------------------------
#  Solve the optmization problem
# ---------------------------------
cd(path_folder_results)
time_to_solve_model = time()                       # Start the timer to solve the Optimization Model
JuMP.optimize!(model)                              # Optimize model
time_to_solve_model = time() - time_to_solve_model # End the timer to build the Optimization Model
println("\nTime to solve the model: $time_to_solve_model sec")
status_model = JuMP.termination_status(model)
println("Termination Status: $status_model \n")
println("--------------------------------------------------------------------------------------------------------------------------------------")
cd(current_path_folder)

#                               ENDS OPTIMIZATION PROCESS 
# ########################################################################################

RBUS::Union{Nothing, DataFrame} = nothing
RGEN::Union{Nothing, DataFrame} = nothing
RCIR::Union{Nothing, DataFrame} = nothing

if status_model == OPTIMAL || status_model == LOCALLY_SOLVED || status_model == ITERATION_LIMIT
    #-------------------------------------------------------------------------------------
    #                             SAVE RESULTS 
    #-------------------------------------------------------------------------------------
    RBUS, RGEN, RCIR = Save_Solution_Model(model, V, θ, P_g, Q_g, P_ik, Q_ik, P_ki, Q_ki, 
    bus_gen_circ_dict_ON, DBUS, DGEN, DGEN_DYN, DCIR, base_MVA, nBUS, nGEN, nCIR, bus_mapping,
    reverse_bus_mapping, current_path_folder, path_folder_results)

    #-------------------------------------------------------------------------------------
    #                             SAVE DUALS 
    #-------------------------------------------------------------------------------------
    Save_Duals_ACOPF_Model(model, V, θ, P_g, Q_g, P_ik, Q_ik, P_ki, Q_ki, eq_const_angle_sw, eq_const_p_balance, eq_const_q_balance, 
    eq_const_p_ik, eq_const_q_ik, eq_const_p_ki, eq_const_q_ki, ineq_const_s_ik, ineq_const_s_ki, 
    ineq_const_diff_ang, base_MVA, current_path_folder, path_folder_results)


    Save_Duals_Dynamic_Model(model, E, δ, Pm, eq_const_Pm, Pe_tf, δ_tf, Δω_tf, δCOI_tf, ΔωCOI_tf, eq_const_δCOI_tf, eq_const_ΔωCOI_tf, eq_const_Pe_tf, eq_const_δ_tf,
    eq_const_Δω_tf, ineq_const_δ_COI_tf, Pe_tpf, δ_tpf, Δω_tpf, δCOI_tpf, ΔωCOI_tpf, eq_const_δCOI_tpf, eq_const_ΔωCOI_tpf, eq_const_Pe_tpf,
    eq_const_δ_tpf, eq_const_Δω_tpf, ineq_const_δ_COI_tpf, base_MVA, current_path_folder, path_folder_results)

    #-------------------------------------------------------------------------------------
    #                          SAVE DYNAMIC RESULTS 
    #-------------------------------------------------------------------------------------
    Manage_Dyn_Results(model, DGEN_DYN, E, δ, Pm, Pe_tf, δ_tf, Δω_tf, δCOI_tf, ΔωCOI_tf, Pe_tpf, δ_tpf,
    Δω_tpf, δCOI_tpf, ΔωCOI_tpf, t_window_total, t_clear_fault, δ_tol, base_MVA, f_syn, ω_syn, current_path_folder, path_folder_results)

else
    JuMP.@warn "Optmization process failed. No feasible solution found."
end
println("--------------------------------------------------------------------------------------------------------------------------------------")

