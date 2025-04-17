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
        mdot=A*p_up/(R*T_up)^0.5*(2*k/(k-1.0)*pr^(2.0/k)*(1-pr^((k-1.0)/k)))^0.5
        #throat temperature
        T_down = T_up * (p_down/p_up)^((k-1.0)/k)
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


#end module
end