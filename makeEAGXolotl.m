% makes xolotl objects, one with
% WT EAG channels and the other with 
% mutant EAG channels 

function [x,x2] =  makeEAGXolotl(self)

p = self.parameters;

x = xolotl;
x.cleanup;
x.closed_loop = false;

x.dt = self.dt;
x.t_end = self.t_end;
T = (self.dt:self.dt:self.t_end);
a = find(T>700,1,'first');
z = find(T>703,1,'first');
V_clamp = 0*(x.dt:x.dt:x.t_end) - 50;
V_clamp(a:z) = 70;

x.addCompartment('C0',-70,p.Ca_in,self.Cm,self.A,self.A,p.phi,self.Ca_out,p.Ca_in,self.tau_Ca,0);
x.V_clamp = V_clamp; % only the first compartment can be clamped 

x.addCompartment('C1',-70,p.Ca_in,self.Cm,self.A,self.A,p.phi,self.Ca_out,p.Ca_in,self.tau_Ca,0);
x.addConductance('C1','aldrich/DmNaV',p.gNa2,self.E_Na);
x.addConductance('C1','Shaker',p.gK1,self.E_K);
x.addConductance('C1','EAGwt',p.gEAG1,self.E_K);
x.addConductance('C1','Leak',p.gLeak,self.E_leak);

x.addCompartment('C2',-70,p.Ca_in,self.Cm,self.A,self.A,p.phi,self.Ca_out,p.Ca_in,self.tau_Ca,0);
x.addConductance('C2','aldrich/DmNaV',p.gNa2,self.E_Na);
x.addConductance('C2','Shaker',p.gK2,self.E_K);
x.addConductance('C2','EAGwt',p.gEAG2,self.E_K);
x.addConductance('C2','hara/Cac',p.gCa2,20);
x.addConductance('C2','Leak',p.gLeak,self.E_leak)

x.addSynapse('Elec','C0','C1',p.g1); 
x.addSynapse('Elec','C0','C2',p.g2); 

x.addSynapse('Elec','C1','C2',p.g12); 
x.addSynapse('Elec','C2','C1',p.g12); 

% make another object for the mutated channels
x2 = xolotl;
x2.closed_loop = false;
x2.dt = self.dt;
x2.t_end = self.t_end;
x2.addCompartment('C0',-70,p.Ca_in,self.Cm,self.A,self.A,p.phi,self.Ca_out,p.Ca_in,self.tau_Ca,0);
x2.V_clamp = V_clamp; % only the first compartment can be clamped 

x2.addCompartment('C1',-70,p.Ca_in,self.Cm,self.A,self.A,p.phi,self.Ca_out,p.Ca_in,self.tau_Ca,0);
x2.addConductance('C1','aldrich/DmNaV',p.gNa2,self.E_Na);
x2.addConductance('C1','Shaker',p.gK1,self.E_K);
x2.addConductance('C1','EAGmut',p.gEAG1,self.E_K);
x2.addConductance('C1','Leak',p.gLeak,self.E_leak);

x2.addCompartment('C2',-70,p.Ca_in,self.Cm,self.A,self.A,p.phi,self.Ca_out,p.Ca_in,self.tau_Ca,0);
x2.addConductance('C2','aldrich/DmNaV',p.gNa2,self.E_Na);
x2.addConductance('C2','Shaker',p.gK2,self.E_K);
x2.addConductance('C2','EAGmut',p.gEAG2,self.E_K);
x2.addConductance('C2','hara/Cac',p.gCa2,20);
x2.addConductance('C2','Leak',p.gLeak,self.E_leak)

x2.addSynapse('Elec','C0','C1',p.g1); 
x2.addSynapse('Elec','C0','C2',p.g1); 

x2.addSynapse('Elec','C1','C2',p.g12); 
x2.addSynapse('Elec','C2','C1',p.g12); 