
using Parameters

@with_kw type synapticTerminalParameters
    
    initial_conditions::Array{Float64,1} = [-50,-50,.05,.1,.1,.1,.1,.1,.1,.1,.1]


    # define default parameters
    cm::Float64 = 10 # nF/mm^2
    A1::Float64 = 0.628e-1 # sq mm
    A2::Float64 = 0.628e-1 # sq mm

    extCa::Float64 = 3000 # in uM
    Ca0::Float64 = 0.05 # uM
    tau_Ca::Float64 = 200 # in milliseconds.
    f::Float64 = 14.96 #microMol/nA

    # electrical coupling
    g1::Float64 = 100;
    g2::Float64 = 100;
    g12::Float64 = 100;

    # conductances
    gNa1::Float64   = 100;
    gNa2::Float64   = 100;
    gCa2::Float64   = 100;
    gK1::Float64    = 100;
    gK2::Float64    = 100;
    g_leak::Float64 = 0;
    

    # reversal potentials
    E_Na::Float64   = 50;
    E_leak::Float64 = -50;
    E_K::Float64    = -80;

end

# define simulation parameters
@with_kw type synTermSimParameters
    dt::Float64 = 50e-3 # ms
    t_end::Int64 = 5000 # ms

    V_drive::Array{Float64,1} = [0.0,0.0] # this, in general, will have to be externally defined. make sure it is the right length (depends on dt, t_end)
end

