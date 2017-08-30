# this file serves as a bridge between MATLAB and Julia
# this file defines various functions that accept MATLAB arguments,
# and run various simulations in Julia and pass the results back to MATLAB


function intSynTerm(neuron_params,reversal_potentials,initial_conditions, time_params, stimulus)
    p = synapticTerminalParameters();
    s = synTermSimParameters();


    # update parameters
    p.cm     = neuron_params[1];
    p.A1     = neuron_params[2];
    p.A2     = neuron_params[3];
    p.extCa  = neuron_params[4]; 
    p.Ca0    = neuron_params[5];
    p.tau_Ca = neuron_params[6];
    p.f      = neuron_params[7];
    p.gNa1   = neuron_params[8];
    p.gNa2   = neuron_params[9];
    p.gCa2   = neuron_params[10];
    p.gK1    = neuron_params[11];
    p.gK2    = neuron_params[12];
    p.g1     = neuron_params[13];
    p.g2     = neuron_params[14];
    p.g12    = neuron_params[15];
    p.g_leak = neuron_params[16];

    p.E_Na = reversal_potentials[1];
    p.E_K  = reversal_potentials[2];
    p.E_leak = reversal_potentials[3];

    p.initial_conditions = initial_conditions;

    s.dt      = Float64(time_params[1]);
    s.t_end   = Int64(time_params[2]);
    s.V_drive = stimulus;


    N = integrate(s,p);
    return N
end