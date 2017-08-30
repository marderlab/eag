% set up the synaptic terminal model -- testing script



s = synapticTerminal;
V = 0*(s.dt:s.dt:s.t_end) -  50;
V(3e4:5e4) = 30;
s.V_drive = V;
N =  s.integrate;