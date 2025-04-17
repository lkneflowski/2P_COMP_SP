using FiniteDifferences
using DifferentialEquations
using Plots
using CoolProp
using ProgressMeter
using Sundials
using DataFrames
using CSV

# include the module valves to get the valves specific functions 
include("valves_single_phase.jl")
# include the state functions
include("state_single_phase.jl")

using .valves_single_phase
using .state_single_phase

# define global variables, series 
global progress_plot = nothing
global time_series = Float64[]
global T_series = Float64[]
global ρ_series = Float64[]
global p_series = Float64[]
global m_series = Float64[]
global y_s_series = Float64[]
global y_d_series = Float64[]
global y_s_saved = Float64[]
global y_d_saved = Float64[]
global saved_y_s = Float64[]
global saved_y_d = Float64[]
global saved_t = Float64[]
global progress_plot
global valve_d_locked_series = Bool[]
global valve_d_stopper_series = Bool[]
global t_valve_series = Float64[]

# working fluid
fluid = "R134a"
# rotational speed
n = 25.0 # Hz
# evap/cond pressure
p_o, p_c = 3e5, 8e5

ω = 2 * pi * n

# superheating in K 
dT_sh = 10 
# calc evap temp 
T_evap = PropsSI("T", "P", p_o, "Q", 1, fluid)
# calc suction temp 
T_0 = PropsSI("T", "P", p_o, "T", T_evap + dT_sh, fluid)

# calc the suction density
ρ_0 = PropsSI("D", "P", p_o, "T", T_0, fluid)
# calc the inlet enthalpy 
h_in = PropsSI("H", "T", T_0, "D", ρ_0, fluid)


#make an initial guess for the density and the temperature --> isentropic compression 
rho_start, T_start = state_single_phase.guess_working_chamber_state(ρ_0, T_0, p_c, fluid)

#adiabatic compressor 
Q_dot = 0.0 

#compressor geometry inputs
bore = 0.030
stroke = 0.033
crank_rad = stroke/2
rod_len = 0.0725

V_disp = pi/4 * bore^2 * stroke
V_clear = V_disp * 0.045
V_tot = V_disp + V_clear

# mass of refrigerant in the cylinder at theta = 0 
m_start = V_clear * rho_start


# valve geometry inputs 
# discharge valve

# youngs modulus of the valve steel
E_d = 2.1e11
# height of the valve
h_valve_d = 0.0003
# valve length 
l_valve_d = 0.0425
af_valve_d = 0.028
#density of the valve steel in kg/m^3
rho_valve_d = 7850
c_w_d = 1.1 
D_valve_d = 0.016
#stopper position 
y_stop_d = 0.00227
b_valve_d = 0.012
I_d = (b_valve_d * h_valve_d^3) / 12
#valve stiffness
k_valve_d = (6 * E_d * I_d) / (af_valve_d^2 * (3 * l_valve_d - af_valve_d))
#mass of the valve 
m_valve_d = 0.0012
# effektive valve mass 
m_eff_d = (1/3) * m_valve_d
#port diameter 
D_port_d = 0.009
#transitional valve lift 
y_tran_d = 0.25 * (D_port_d^2 / D_valve_d)
#port area
A_port_d = pi * (D_port_d^2)/4
#valve area 
A_valve_d = pi *(D_valve_d^2)/4


# suction valve 
# youngs modulus of the valve steel 
E_s = 2.1e11
# height of the valve
h_valve_s = 0.0003
# valve length 
l_valve_s = 0.044
af_valve_s = 0.0295
#density of the valve steel in kg/m^3
rho_valve_s = 7850
c_w_s = 1.1 
D_valve_s = 0.016
#stopper position 
y_stop_s = 0.0004065

b_valve_s = 0.012
I_s = (b_valve_s * h_valve_s^3) / 12
#valve stiffness
k_valve_s = (6 * E_s * I_s) / (af_valve_s^2 * (3 * l_valve_s - af_valve_s))
#mass of the valve 
m_valve_s = 0.0012
# effektive valve mass 
m_eff_s = (1/3) * m_valve_s
#port diameter 
D_port_s = 0.012
#transitional valve lift 
y_tran_s = 0.25 * (D_port_s^2 / D_valve_s)
#port area
A_port_s = pi * (D_port_s^2)/4
#valve area 
A_valve_s = pi *(D_valve_s^2)/4



#initial values for the valves 
#suction valve initial values
y_s_start = 0.0
v_s_start = 0.0 
#discharge valve initial values 
y_d_start = 0.0
v_d_start = 0.0 



function s_function(r, l, ω)
    s(t) = r * (1 + l/r - cos(ω * t ) - (l/r) * sqrt(1 - (r/l)^2 * sin(ω * t)^2))
    return s
end

function V_function(s_fun, bore, V_clear)
    V(t) = s_fun(t) * pi/4 * bore^2 + V_clear
    return V
end


dVdt(t) = central_fdm(5, 1)(V_t, t)
∂p∂T(T, rho, fluid) = central_fdm(5, 1)(T->PropsSI("P", "T", T, "D", rho, fluid), T)
c_v(T, rho, fluid) = PropsSI("CVMASS", "T", T, "D", rho, fluid)
h(T, rho, fluid) = PropsSI("H", "T", T, "D", rho, fluid)


s_t = s_function(crank_rad, rod_len, ω)
V_t = V_function(s_t, bore, V_clear)
#ṁ_p = mass_flow(p_o, p_c, ρ_0, fluid)
∂p∂T_s(u) = ∂p∂T(u[2], u[1], fluid)
c_vu(u) = c_v(u[2], u[1], fluid)
h_u(u) = h(u[2], u[1], fluid)
#Ḣ_u = H_flow(h_in, fluid)

Q_dot = 0 

mutable struct CompressorParams
    suction_valve::NTuple{7, Float64}
    discharge_valve::NTuple{7, Float64}
    common::NTuple{6, Any}
    valve_s_locked::Bool
    valve_d_locked::Bool
    valve_s_stopper:: Bool
    valve_d_stopper:: Bool 
end

# Initialize parameters
p = CompressorParams(
    (c_w_s, A_valve_s, A_port_s, k_valve_s, m_eff_s, y_tran_s, y_stop_s),
    (c_w_d, A_valve_d, A_port_d, k_valve_d, m_eff_d, y_tran_d, y_stop_d),
    (ρ_0, T_0, fluid, ω, p_c, Q_dot),
    false,
    false, 
    false, 
    false
)


function compressor!(du, u, p, t)

    # Entpacken der Parameter
    # Parameter für das Saugventil
    c_w_s, A_valve_s, A_port_s, k_valve_s, m_eff_s, y_tran_s, y_stop_s = p.suction_valve
    c_w_d, A_valve_d, A_port_d, k_valve_d, m_eff_d, y_tran_d, y_stop_d = p.discharge_valve
    ρ_0, T_0, fluid, ω, p_c, Q_dot = p.common
    valve_s_locked = p.valve_s_locked
    valve_d_locked = p.valve_d_locked
    valve_s_stopper = p.valve_s_stopper
    valve_d_stopper = p.valve_d_stopper

    

    #u[1] = T 
    #u[2] = ρ
    #u[3] = m 
    y_s  = u[4]
    v_s = u[5]
    y_d = u[6]
    v_d = u[7]

    #VALVES 
    #upstreamstate infront of the suction valve (suction conditions)
    State_up_s = state_single_phase.calc_state(ρ_0, T_0, fluid)
    # calc the state in the working chamber from the solution of the thermodynamic model
    State_down_s = state_single_phase.calc_state(u[1], u[2], fluid)
    # Berechnung der Strömungsgeschwindigkeit durch das Saugventil
    w_t_s, mdot_s = valves_single_phase.flow_velocity(u[4], y_tran_s, State_up_s, State_down_s, D_valve_s, A_port_s)
    

    State_up_d = State_down_s
    State_down_d = Dict("p" => p_c)

    w_t_d, mdot_d = valves_single_phase.flow_velocity(u[6], y_tran_d, State_up_d, State_down_d, D_valve_d, A_port_d)


    #WORKING CHAMBER 
    #drho/dt
    du[1] = 1/V_t(t) * (-u[1] * dVdt(t) + du[3])
    #dT/dt
    du[2] = (-u[2] * (∂p∂T_s(u)) * (dVdt(t) - 1/u[1] * du[3]) - h_u(u) * du[3] + Q_dot + (mdot_s * h(T_0, ρ_0, fluid) - mdot_d * h_u(u)))  / (u[3] * c_vu(u))
    #dm/dt
    du[3] = mdot_s - mdot_d


    
    if y_s > y_tran_s # flux driven region
        du[4] = v_s  # dy/dt = v
        du[5] = (1/m_eff_s) * ( (1/2) * c_w_s * State_up_s["rho"] * w_t_s^2 * A_valve_s + State_up_s["rho"] * (w_t_s - v_s)^2 * A_port_s - k_valve_s * u[4])
        
    else          # pressure driven region 
        du[4] = v_s  # dy/dt = v
        du[5] = (1/m_eff_s) * ( (1/2) * c_w_s * State_up_s["rho"] * w_t_s^2 * A_valve_s + (State_up_s["p"] - State_down_s["p"])*A_valve_s - k_valve_s * u[4])
    end 
    


    
    if y_d > y_tran_d # flux driven region
        du[6] = v_d  # dy/dt = v
        du[7] = (1/m_eff_d) * ( (1/2) * c_w_d * State_up_d["rho"] * w_t_d^2 * A_valve_d + State_up_d["rho"] * (w_t_d - v_d)^2 * A_port_d - k_valve_d * u[6])
        
    else          # pressure driven region 
        du[6] = v_d  # dy/dt = v
        du[7] = (1/m_eff_d) * ( (1/2) * c_w_d * State_up_d["rho"] * w_t_d^2 * A_valve_d + (State_up_d["p"] - State_down_d["p"])*A_valve_d - k_valve_d * u[6])
    end 

    if p.valve_s_locked
        u[8] = 0.0 
    
    elseif p.valve_s_stopper
        u[8] = y_stop_s
    else 
        u[8] = u[4]
    end


    if p.valve_d_locked
        u[9] = 0.0 

    elseif p.valve_d_stopper
        u[9] = y_stop_d
    else
        u[9] = u[6]
    end

end



function dynamic_plot_callback(u, t, integrator)
    # Declare global variables
    global progress_plot, time_series, T_series, ρ_series, p_series, m_series, y_s_series, y_d_series

    # Current time and state from function arguments
    current_time = t
    current_state = u
    

    # Extract specific state variables
    T = current_state[2]  # Temperature
    ρ = current_state[1]  # Density
    m = current_state[3]
    y_s = current_state[4]
    y_d = current_state[6]

    # Calculate pressure based on temperature and density
    fluid = "R134a"
    p = PropsSI("P", "D", ρ, "T", T, fluid)

    # Append values to time series
    push!(time_series, current_time)
    push!(T_series, T)
    push!(ρ_series, ρ)
    push!(p_series, p)
    push!(m_series, m)
    push!(y_s_series, y_s)
    push!(y_d_series, y_d)

    # Initialize the plot only once
    if progress_plot === nothing
        progress_plot = plot(
            layout=(6, 1),
            size=(1000, 900)
        )
        # Create subplots without legends (they update automatically)
        plot!(progress_plot[1], time_series, T_series, label="", xlabel="Zeit [s]", ylabel="Temperatur [K]")
        plot!(progress_plot[2], time_series, ρ_series, label="", xlabel="Zeit [s]", ylabel="Dichte [kg/m³]")
        plot!(progress_plot[3], time_series, p_series, label="", xlabel="Zeit [s]", ylabel="Druck [Pa]")
        plot!(progress_plot[4], time_series, m_series, label="", xlabel="Zeit [s]", ylabel="Masse m")
        plot!(progress_plot[5], time_series, y_s_series, label="", xlabel="Zeit [s]", ylabel="Ventilhub Saugventil")
        plot!(progress_plot[6], time_series, m_series, label="", xlabel="Zeit [s]", ylabel="Ventilhub Druckventil")
    else
        # Update the subplots with new data
        plot!(progress_plot[1], time_series, T_series, label="", overwrite = true)
        plot!(progress_plot[2], time_series, ρ_series, label="",overwrite = true)
        plot!(progress_plot[3], time_series, p_series, label="", overwrite = true)
        plot!(progress_plot[4], time_series, m_series, label="", overwrite = true)
        plot!(progress_plot[5], time_series, y_s_series, label="", overwrite = true, ylim=(-0.005, 0.002))
        plot!(progress_plot[6], time_series, y_d_series, label="", overwrite = true)
    end

    # Display the updated plot
    display(progress_plot)
end

#SUCTIONVALVE
#__________________________________________________________________________________________

# Callback 1 --> Discrete when suction valve lift ist smaller tahn 0.0 
function condition_y_zero_s(u, t, integrator)
    return u[4] < 0  # Bedingung erfüllt, wenn y_s < 0
end

function affect_y_zero_s!(integrator)
    integrator.p.valve_s_locked = true 
    integrator.u[4] = 0.0  # Setze den Ventilhub auf 0
    integrator.u[5] = 0.0  # Geschwindigkeit auf 0 setzen
end

# Callback 2 --> Discrete unlock the valve, when valve lift gets higher than 0.0 
function condition_y_unlock_s(u, t, integrator)
    return integrator.p.valve_s_locked && u[4] > 0
end

function affect_y_unlock_s!(integrator)
    integrator.p.valve_s_locked = false 
end


#Callback 3 --> Discrete set the valve to stopper position if the valve lift is greater than the stopper 
function condition_y_stop_s(u, t, integrator)
    return u[4] > y_stop_s  # Bedingung erfüllt, wenn y_s > y_stop_s
end

function affect_y_stop_s!(integrator)
    integrator.p.valve_s_stopper = true 
    integrator.u[4] = y_stop_s  # Setze den Ventilhub auf den maximalen Wert
    integrator.u[5] = 0.0       # Geschwindigkeit auf 0 setzen
end

#Callback 4 --> Discrete unlock the valve, when the valve lift is again smaller than y_stop_s
function condition_y_unlock_s_stopper(u, t, integrator)
    return integrator.p.valve_s_stopper && u[4] < y_stop_s
end

function affect_y_unlock_s_stopper!(integrator)
    integrator.p.valve_s_stopper = false 
end


# Callback 5 --> Continous Callback to detect the transfer (to be checked)
function condition_y_zero_s_transfer(u, t, integrator)
    return u[4] 
    
end 

function affect_y_zero_s_transfer!(integrator)
    integrator.p.valve_s_locked = true 
    integrator.u[4] = 0.0  # Setze den Ventilhub auf 0
    #integrator.u[5] = 0.0  # Geschwindigkeit auf 0 setzen
end 

#Callback6 --> Continous Callback to detect the transfer (to be checked)
function condition_y_stop_s_transfer(u, t, integrator)
    return u[4] - y_stop_s# Bedingung erfüllt, wenn y_s > y_stop_s
end

function affect_y_stop_s_transfer!(integrator)
    integrator.p.valve_s_stopper = true 
    integrator.u[4] = y_stop_s  # Setze den Ventilhub auf den maximalen Wert
    #integrator.u[5] = 0.0       # Geschwindigkeit auf 0 setzen
end
#__________________________________________________________________________________________


#DISCHARGE VALVE
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Callback 1 --> Discrete when discharge valve lift ist smaller than 0.0 
function condition_y_zero_d(u, t, integrator)
    return u[6] < 0  # Bedingung erfüllt, wenn y_d < 0
end

function affect_y_zero_d!(integrator)
    integrator.p.valve_d_locked = true
    integrator.u[6] = 0.0  # Setze den Ventilhub auf 0
    integrator.u[7] = 0.0  # Geschwindigkeit auf 0 setzen
end


# Callback 2 --> Discrete unlock the valve, when valve lift gets higher than 0.0
function condition_y_unlock_d(u, t, integrator)
    return integrator.p.valve_d_locked && u[6] > 0
end

function affect_y_unlock_d!(integrator)
    integrator.p.valve_d_locked = false
end



#Callback 3 --> Discrete set the valve to stopper position if the valve lift is greater than the stopper 
function condition_y_stop_d(u, t, integrator)
    return u[6] > y_stop_d  # Bedingung erfüllt, wenn y_d > y_stop_d
end

function affect_y_stop_d!(integrator)
    integrator.p.valve_d_stopper = true 
    integrator.u[6] = y_stop_d  # Setze den Ventilhub auf den maximalen Wert
    integrator.u[7] = 0.0       # Geschwindigkeit auf 0 setzen
end

#Callback 4 --> Discrete unlock the valve, when the valve lift is again smaller than y_stop_s
function condition_y_unlock_d_stopper(u, t, integrator)
    return integrator.p.valve_d_stopper && u[6] < y_stop_d
end

function affect_y_unlock_d_stopper!(integrator)
    integrator.p.valve_d_stopper = false 
end



#Callback 5 --> Continous to be checked
function condition_y_zero_d_transfer(u, t, integrator)
    return u[6] 
    
end 

function affect_y_zero_d_transfer!(integrator)
    integrator.p.valve_d_locked = true 
    integrator.u[6] = 0.0  # Setze den Ventilhub auf 0
    #integrator.u[7] = 0.0  # Geschwindigkeit auf 0 setzen
end 

function condition_y_stop_d_transfer(u, t, integrator)
    return u[6] - y_stop_d # Bedingung erfüllt, wenn y_d > y_stop_d
end

function affect_y_stop_d_transfer!(integrator)
    integrator.p.valve_d_stopper = true 
    integrator.u[6] = y_stop_d  # Setze den Ventilhub auf den maximalen Wert
    #integrator.u[7] = 0.0       # Geschwindigkeit auf 0 setzen
end



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function log_valve_states(u, t, integrator)
    push!(valve_d_locked_series, integrator.p.valve_d_locked)
    push!(valve_d_stopper_series, integrator.p.valve_d_stopper)
    push!(t_valve_series, t)
end


# CallbackSet erstellen
function create_callbacks()

    #Callback 1 - lock the valve discrete for y < 0.0 
    callback_y_zero_s = DiscreteCallback(condition_y_zero_s, affect_y_zero_s!, save_positions = (true, true))
    #Callback 2 - unlock the valve when y > 0.0 
    callback_unlock_s = DiscreteCallback(condition_y_unlock_s, affect_y_unlock_s!, save_positions = (true, true) )
    #Callback 3 - stop the valve when y > y_stop_s
    callback_y_s = DiscreteCallback(condition_y_stop_s, affect_y_stop_s!, save_positions = (true, true))


    callback_y_zero_s_transfer = ContinuousCallback(condition_y_zero_s_transfer, affect_y_zero_s_transfer!, save_positions = (true, true))
    #Callback 4 - unlock the valve when y < y_stop_s again (downcrossing from y_stop_s)
    #upcrossing valve lift --> lock the valve continously
    callback_y_s_transfer = ContinuousCallback(condition_y_stop_s_transfer, affect_y_stop_s_transfer!, nothing, save_positions = (true, true) )
    callback_unlock_s_stopper = DiscreteCallback(condition_y_unlock_s_stopper, affect_y_unlock_s_stopper!,save_positions = (true, true) )
    




    callback_y_zero_d = DiscreteCallback(condition_y_zero_d, affect_y_zero_d!, save_positions = (true, true))
    callback_unlock_d = DiscreteCallback(condition_y_unlock_d, affect_y_unlock_d!, save_positions = (true, true) )
    callback_y_d = DiscreteCallback(condition_y_stop_d, affect_y_stop_d!, save_positions = (true, true))
    callback_y_zero_d_transfer = ContinuousCallback(condition_y_zero_d_transfer, affect_y_zero_d_transfer!, save_positions = (true, true))
    callback_y_d_transfer = ContinuousCallback(condition_y_stop_d_transfer, affect_y_stop_d_transfer!, nothing, save_positions = (true, true))
    callback_unlock_d_stopper = DiscreteCallback(condition_y_unlock_d_stopper, affect_y_unlock_d_stopper!, save_positions = (true, true) )

    #log_valves_cb = FunctionCallingCallback(log_valve_states)

    return CallbackSet(callback_y_zero_s, callback_unlock_s, callback_y_s, callback_unlock_s_stopper, callback_y_s_transfer, callback_y_s_transfer,      
                       callback_y_zero_d, callback_unlock_d, callback_y_d, callback_y_zero_d_transfer, callback_unlock_d_stopper, callback_y_d_transfer)
end

global last_progress = 0.0
global t_end = 2*pi/ω

function progress_callback(u, t, integrator)
    progress = t / t_end * 100
    if progress - last_progress >= 10  # Zeige Fortschritt in 10%-Schritten
        println("Progress: $(round(progress, digits=1))%")
        global last_progress = progress
    end
end

progress_cb = FunctionCallingCallback(progress_callback)

callbacks = create_callbacks()

plot_callback = FunctionCallingCallback(dynamic_plot_callback)

#plot_callback = SavingCallback((u, t, integrator) -> dynamic_plot_callback(integrator),save_positions=(false, false))

# Kombiniere mit anderen Callbacks
combined_callbacks = CallbackSet(callbacks, progress_cb)



u₀ = [rho_start, T_start, m_start, y_s_start, v_s_start, y_d_start, v_d_start, y_s_start, y_d_start]


# Zeitspanne für φ ∈ [0, π/4]
tspan = (0.0, 2π / (1 * ω))




prob = ODEProblem(compressor!, u₀, tspan, p)



#sol = solve(prob, Tsit5(), callback = combined_callbacks, reltol=1e-7, abstol=1e-7, dtmax = 1e-3)
sol = solve(prob, BS3(), callback = combined_callbacks, reltol=1e-4, abstol=1e-4)

plot(sol.t, sol.u[:, 1], label="ρ(t)", xlabel="Zeit", ylabel="ρ")

sol_matrix = hcat(sol.u...)'
sol_rho = sol_matrix[:, 1]
sol_T = sol_matrix[:, 2]
sol_m = sol_matrix[:,3]
sol_y_s = sol_matrix[:,4]
sol_y_d = sol_matrix[:,6]
sol_y_s_dummy = sol_matrix[:,8]
sol_y_d_dummy = sol_matrix[:,9]
dt = diff(sol.t)
#φ_mid = φ[1:end-1] .+ diff(φ) ./ 2  # Kurbelwinkelwerte in der Mitte der Schritte

sol_p = PropsSI.("P", "T", sol_T, "D", sol_rho, fluid)
sol_s = PropsSI.("S", "T", sol_T, "D", sol_rho, fluid)

φ = sol.t .* ω  # φ(t) = ω * t
φ_valve = ω .* t_valve_series

# DataFrame erstellen
data = DataFrame(
    theta = φ,
    rho = sol_rho,
    T = sol_T,
    m = sol_m,
    y_s = sol_y_s,
    y_d = sol_y_d,
    y_s_dummy = sol_y_s_dummy,
    y_d_dummy = sol_y_d_dummy,
    p = sol_p,
    s = sol_s,
)


file_path = "C:\\Users\\Leonard Kneflowski\\PycharmProjects\\2P-Compresion\\simulation_results_main.csv"


# Datei speichern
#CSV.write(file_path, data)



# Plot des Drucks über den Kurbelwinkel
plot(
plot(φ, sol_p*1e-5, label="Druck p(φ)", xlabel="Kurbelwinkel φ [rad]", ylabel="Druck p [Pa]", legend=:topleft),
plot(φ, sol_m, xlabel="Kurbelwinkel φ [rad]", ylabel="Masse m", legend=:topleft),
plot(φ, sol_y_s, xlabel="Kurbelwinkel φ [rad]", ylabel="Ventilhub S", legend=:topright),
plot(φ, sol_y_s_dummy, xlabel="Kurbelwinkel φ [rad]", ylabel="Ventilhub S", legend=:topright),
plot(φ, sol_y_d, xlabel="Kurbelwinkel φ [rad]", ylabel="Ventilhub D", legend=:topright),
plot(φ, sol_y_d_dummy, xlabel="Kurbelwinkel φ [rad]", ylabel="Ventilhub D", legend=:topright),
#plot(φ_valve, valve_d_locked_signal, label="D-Ventil locked", xlabel="Kurbelwinkel φ [rad]", ylabel="Zustand", ylim=(-0.1, 1.1)),
#plot(φ_valve, valve_d_stopper_signal, label="D-Ventil Stopper", xlabel="Kurbelwinkel φ [rad]", ylabel="Zustand", ylim=(-0.1, 1.1)),
#plot(φ, dt, xlabel="Kurbelwinkel φ [rad]", ylabel="Schrittweite Δt", label="Solver Schrittweite", legend=:topright),
layout=(3, 2),  # 4 Zeilen, 2 Spalten
size=(1200, 1000),
)


