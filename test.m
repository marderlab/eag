% set up the synaptic terminal model -- testing script

s = synapticTerminal;
s.t_end = 500;
time = s.dt:s.dt:s.t_end;

V = 0*time -  50;
a = find(time>100,1,'first');
z = find(time>103,1,'first'); % 3 ms wide AP
V(a:z) = 50;


s.V_drive = V;
s.manipulate;



%% for voltage clamp
s = synapticTerminal;
s.t_end = 1e3;
time = s.dt:s.dt:s.t_end;

V = 0*time -  51;
V(200:end-201) = 10;

s.V_drive = V;
s.manipulate;
