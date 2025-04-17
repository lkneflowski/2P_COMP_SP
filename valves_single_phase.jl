#file that contains the functions of the valve model 

module valves_single_phase

export A_flow, flow_velocity, isentropic_nozzle, valve_dynamics


# function to calculate the flow area A (implmented from PDSim)
function A_flow(y, y_tran, D_valve, A_port)
    if y === nothing 
        #println("y = None")
    end
    if y >= y_tran
        return A_port
    else 
        return pi * y * D_valve
    end
end
 
#function to calculate the valve flow velocity 
function flow_velocity(y, y_tran, State_up, State_down, D_valve, A_port)
    
    A = A_flow(y, y_tran, D_valve, A_port)

   # println(A,"A")

    
    if A > 0
        w_curr, m_dot_curr = isentropic_nozzle(A, State_up, State_down)

        if y > y_tran
            return w_curr, m_dot_curr
        else
            return (y/y_tran) * w_curr, m_dot_curr
        end
        #println(m_dot_curr)
    else
        return 0.0, 0.0
    end
end




# function to calculate the isentropic nozzle flow --> returns the flow velocity
function isentropic_nozzle(A, State_up, State_down)
    
    #println(A, "A")
    cp = State_up["cp"]  #J/kg/K
    R = 8.314462618/State_up["M"]      #J/mol K * (mol/kg) = J/K*kg
    cv = cp - R

    k = cp/cv 
    p_up = State_up["p"]
    T_up = State_up["T"]
    p_down = State_down["p"]
    #speed of sound 
    c = (k*R*T_up)^(0.5)
    #upstream density 

    rho_up = p_up /(R*T_up)   # Pa /( J/K*kg * K )  --> kg/m^3
    #println(rho_up,"rho_up")
    pr = p_down/p_up
    pr_crit = (1+ (k-1)/2.0)^(k/(1-k))

   # println(cp,"cp")
    #println(R, "R")
   # print(cv,"cv")
    #println(k,"k")
   # println(p_up,"p_up")
   # println(p_down,"p_down")
    #println(T_up, "T_up")
    #println(pr, "pressure_ratio")
    #println((pr^(2/k) - pr^((k+1)/k)), "Ausdruck")


    if pr > 1  
        # no flow through the valves if downstream pressure is higher than upstream pressure!
        w = 0 
        mdot = 0 
        return w, mdot 
    end

    #wenn das kritische Druckverhältnis überschritten wird 
    if pr > pr_crit
        #mass flow rate is not choked
        #mdot = A * sqrt(p_up * rho_up) * sqrt(((2*k/(k-1))*(pr^(2/k) - pr^((k+1)/k))))
        mdot = (A * p_up /(sqrt(R * T_up))) * (2*k/(k-1.0))*pr^(2.0/k) * sqrt((1-pr^((k-1)/k)))
        #throat temperature
        T_down = T_up * (p_down/p_up)^((k-1)/k)
        #throat density 
        rho_down = p_down/(R*T_down)
        #Velocity at throat 
        w = mdot/(rho_down*A)
        # Mach number 
        Ma = w/c
    else
        #choked mass flow rate 
        mdot = A*rho_up*(k*R*T_up)^0.5*(1+(k-1)/2)^((1+k)/(2*(1-k)))
        #velocity at throat 
        w = c 
        # Mach number 
        Ma = 1
    end

    return w, mdot


    #println("w", w, "m_dot", mdot)
end


function valve_dynamics(y, v, params, State_up, State_down, w_t)
    c_w, A_valve, A_port, k_valve, m_eff, y_stop, y_tran = params

    # Ventilbeschränkungen
    if y > y_stop
        y = y_stop
    elseif y < 0.0
        y = 0.0
    elseif y < 1e-15
        y = 0.0
    end

    # Ventildynamik
    if y > y_tran  # Strömungsgetriebener Bereich
        dy = v
        dv = (1/m_eff) * (
            0.5 * c_w * State_up["rho"] * w_t^2 * A_valve +
            State_up["rho"] * (w_t - v)^2 * A_port -
            k_valve * y
        )
    else  # Druckgetriebener Bereich
        dy = v
        dv = (1/m_eff) * (
            0.5 * c_w * State_up["rho"] * w_t^2 * A_valve +
            (State_up["p"] - State_down["p"]) * A_valve -
            k_valve * y
        )
    end

    return dy, dv, y
end


#end module
end