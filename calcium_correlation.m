
close all

% demonstrates correlation between regular bursting and <[Ca]>

% instantiate a xolotl model
% >> xolotl object with 
% ---------------------
% + AB  
%   > ACurrent (g=500, E=-80)
%   > CaS (g=60, E=49)
%   > CaT (g=25, E=49)
%   > HCurrent (g=0.1, E=-20)
%   > KCa (g=50, E=-80)
%   > Kd (g=1000, E=-80)
%   > Leak (g=1e-05, E=-55)
%   > NaV (g=1000, E=50)
% ---------------------
x = xolotl.examples.BurstingNeuron('prefix','prinz');
x.AB.add('Leak','gbar',1e-5);
x.t_end = 5e3;

% extract voltage, calcium, mRNA, and current
% `V` Voltage trace of every compartment. A matrix of size (nsteps, n_comps)
% `Ca` Calcium concentration in every cell and the corresponding `E_Ca` (reversal potential of Calcium). A matrix of size (nsteps, n_comps)
% `M` a matrix representing every dimension of every mechanism in the tree. This matrix has size (nsteps, NC), where NC depends on the precise controllers used, and is automatiCally determined.
% `I` the currents of every ion channel type in the model. This is a matrix of size (nsteps, n_cond)
[V,Ca,M,I] = x.integrate;

% set calcium target to average calcium concentration
x.AB.Ca_target = mean(Ca(:,1));

% calculate total calcium current by summing CaS and CaT
iCa = sum(I(:,2:3),2);

% instantiate an exponential low pass filter
K = exp(linspace(0,10000,10000)/-3000);

% plot filtered, normalized calcium current against calcium concentration
figure('outerposition',[300 300 900 600],'PaperUnits','points','PaperSize',[1200 600]); hold on
plot(filter(K,sum(K),iCa)); hold on
plot(Ca(:,1))

drawnow

figlib.pretty('PlotLineWidth',1.5,'LineWidth',1.5)