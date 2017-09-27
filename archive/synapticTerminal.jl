# 3-compartment synaptic terminal, 
# based loosely on a Hodgkin-Huxley neuron 


function integrate(simp::synTermSimParameters, param::synapticTerminalParameters)
    # unpack parameters
    dt::Float64 = simp.dt;
    t_end::Float64 = simp.t_end;
    V_drive::Array{Float64, 1} = simp.V_drive;

    f::Float64 = param.f;
    A1::Float64 = param.A1;
    A2::Float64 = param.A2;
    extCa::Float64 = param.extCa;

    # conductances
    gNa1::Float64 = param.gNa1;
    gK1::Float64  = param.gK1;
    gNa2::Float64 = param.gNa2;
    gK2::Float64  = param.gK2;
    gCa2::Float64 = param.gCa2;
    g_leak::Float64 = param.g_leak;

    g1::Float64 = param.g1;
    g2::Float64 = param.g2;
    g12::Float64 = param.g12;

    # reversal potentials
    E_Na::Float64   = param.E_Na;
    E_leak::Float64 = param.E_leak;
    E_K::Float64    = param.E_K;

    Ca0::Float64 = param.Ca0;
    cm::Float64 = param.cm;


    # preallocate some arrays + initial_conditions
    nsteps::Int64 = floor(t_end/dt) + 1; 
    V1   = Array{Float64, 1}(nsteps); 
    V2   = Array{Float64, 1}(nsteps);
    Ca   = Array{Float64, 1}(nsteps);  

    mNa1::Float64 = rand();
    mK1 ::Float64 = rand();
    hNa1::Float64 = rand();
    mNa2::Float64 = rand();
    mK2 ::Float64 = rand();
    mCa2::Float64 = rand();
    hNa2::Float64 = rand();
    hCa2::Float64 = rand();



    # temporary variables 
    g = Array{Float64, 1}(7);
    zinf = Array{Float64, 1}(11);
    tz = Array{Float64, 1}(11);
    exp_dt_by_tau_Ca::Float64 = exp(-dt/param.tau_Ca);

    # (in)activation variables for the 7 channels (placeholders)
    h = ones(7);
    m = ones(7);
    m_inf = ones(7);
    h_inf = ones(7);
    tau_m = ones(7);
    tau_h = ones(7);
    
    R_by_zF = 500.0*(8.6174e-5);
    T = 10 + 273.15;
    RT_by_zF::Float64 = R_by_zF*T;
    fA1 = f*A1;
    fA2 = f*A2;

    for int_step = 2:nsteps

        V_ext   = V_drive[int_step-1];
        V_prev1 = V1[int_step-1];
        V_prev2 = V2[int_step-1];
        Ca_prev = Ca[int_step-1];

        # find m_inf, h_inf
        m_infNa1 = mNainf(V_prev1);
        h_infNa1 = hNainf(V_prev1);

        m_infK1 = mKdinf(V_prev1);

        # integrate m 
        mNa1 = m_infNa1 + (mNa1 - m_infNa1)*exp(-(dt/taumNa(V_prev1)));
        mK1  = m_infK1 + (mK1 - m_infK1)*exp(-(dt/taumNa(V_prev1)));

        # integrate h
        hNa1 = m_infNa1 + (hNa1 - m_infNa1)*exp(-(dt/tauhNa(V_prev1)));

        # compute currents
        I1 = g1*(V_ext - V_prev1)*1e-6; # in nA

        # compute effective conductances 
        g[1] = mNa1*mNa1*mNa1*hNa1; wtf?? g_bars are missing 
        g[2] = mK1*mK1*mK1*mK1;  
   

        sigma_g = sum(g);

              
        # compute V_inf, tau_V using g <---- m(t), h(t) 
        V_inf1::Float64 = (g[1]*E_Na + g[2]*E_K + g_leak*E_leak + I1/A1)/sigma_g;
        tau_v1::Float64 = cm/(sigma_g*10);

        # integrate V
        @fastmath V1[int_step] = V_inf1 + (V_prev1 - V_inf1)*exp(-(dt/tau_v1)); # mV


    end 

    V1 = V1[2:end]::Array{Float64,1};
    V2 = V2[2:end]::Array{Float64,1};
    Ca = Ca[2:end]::Array{Float64,1};
    n = Array{Float64,2}(nsteps-1,11);
    n[:,1] = V1;
    n[:,2] = V2;
    n[:,3] = Ca;
    return n
end 





