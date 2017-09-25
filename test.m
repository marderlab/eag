

% set up the synaptic terminal model using xolotl 


term_area = 1.25e-5; % mm^2
f = 14.96;
Ca_out = 1500; 
Ca_in = 0.036; 
tau_Ca = 85;
Cm = 10;

x = xolotl;
x.dt = 50e-3;
x.t_end = 1e3;
V_clamp = 0*(x.dt:x.dt:x.t_end) - 50;
V_clamp(5e3:5.1e3) = 70;
x.addCompartment('C1',-70,Ca_in,Cm,term_area,f,Ca_out,Ca_in,tau_Ca);
x.V_clamp = V_clamp; % only the first compartment can be clamped 

x.addCompartment('C2',-70,Ca_in,Cm,term_area,f,Ca_out,Ca_in,tau_Ca);
x.addConductance('C2','prinz/NaV',200,50);
x.addConductance('C2','prinz/Kd',200,-80);

x.addCompartment('C3',-70,Ca_in,Cm,term_area,f,Ca_out,Ca_in,tau_Ca);
x.addConductance('C3','prinz/NaV',200,50);
x.addConductance('C3','prinz/Kd',200,-80);


x.addSynapse('Elec','C1','C2',.05); % one-way electrical syanpse from 1->2
x.addSynapse('Elec','C1','C3',.05); % one-way electrical syanpse from 1->3

x.addSynapse('Elec','C2','C3',.05); 
x.addSynapse('Elec','C3','C2',.05); 

x.compile;

x.manipulate;

