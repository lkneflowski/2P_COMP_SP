module state_single_phase
using CoolProp

export calc_state, calc_downstream_pressure, guess_working_chamber_state


# calculates the initial state for a given rho and T 
function calc_state(rho, T, fluid)

    state = Dict()

    state["T"] = T  #K
    state["rho"] = rho  #kg/m^3
    
    state["p"]= PropsSI("P", "D", state["rho"], "T", state["T"], fluid)   #Pa
    state["h"] = PropsSI("H", "D", state["rho"], "T", state["T"], fluid)  #J/kg
    state["cp"] = PropsSI("CP0MASS", "D", state["rho"], "T", state["T"], fluid) # J/kg/K
    state["s"] = PropsSI("S", "D", state["rho"], "T", state["T"], fluid)   
    #state["cv"] = PropsSI("CVMASS", "D", state["rho"], "T", state["T"], fluid)
    state["M"] = PropsSI("MOLARMASS", "D", state["rho"], "T", state["T"], fluid) #kg/mol

    return state
end

function calc_downstream_pressure(p_cond)
    downstream_pressure = Dict()
    # downstream pressure is the condensing pressure 
    downstream_pressure["p"] = p_cond

    return downstream_pressure
end

function guess_working_chamber_state(rho_suc, T_suc, p_cond, fluid)
    #from rho and T calc entropy s 
    s = PropsSI("S", "T", T_suc, "D", rho_suc, fluid)
    # calc density at entropy s and condensing pressure 
    rho_start = PropsSI("D", "S", s, "P", p_cond, fluid)
    # calc temp
    T_start = PropsSI("T", "S", s, "P", p_cond, fluid)

    return rho_start, T_start

end



end 