% set up the synaptic terminal model -- testing script

s = synapticTerminal;
s.t_end = 1000;
time = s.dt:s.dt:s.t_end;

V = 0*time -  50;
a = find(time>200,1,'first');
z = find(time>300,1,'first');
V(a:z) = 10;
s.V_drive = V;
s.manipulate;

