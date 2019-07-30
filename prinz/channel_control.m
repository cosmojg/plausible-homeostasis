% reproduces Fig 3c in Tim O'Leary's paper
% controltype :
% ==1         : all concentration integral
% ==2         : all current integral
% >=3         : random mix
function x = channel_control(controltype)

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
% x = xolotl.examples.BurstingNeuron('prefix','liu');
x = xolotl.examples.BurstingNeuron('prefix','prinz');
x.sim_dt = .1;
x.dt = .1;

% extract voltage, calcium, mRNA, and current
% `V` Voltage trace of every compartment. A matrix of size (nsteps, n_comps)
% `Ca` Calcium concentration in every cell and the corresponding `E_Ca` (reversal potential of Calcium). A matrix of size (nsteps, n_comps)
% `M` a matrix representing every dimension of every mechanism in the tree. This matrix has size (nsteps, NC), where NC depends on the precise controllers used, and is automatiCally determined.
% `I` the currents of every ion channel type in the model. This is a matrix of size (nsteps, n_cond)
x.t_end = 1e5;
[~,~,M,~] = x.integrate;

% set calcium target to average calcium concentration
x.AB.Ca_target = x.AB.Ca_average;

% set channel controllers
controllers = {'IntegralController', 'IntegralCurrentController'};
if controltype == 1
    x.AB.ACurrent.add(strjoin(controllers(1)));
    x.AB.CaS.add(strjoin(controllers(1)));
    x.AB.CaT.add(strjoin(controllers(1)));
    x.AB.HCurrent.add(strjoin(controllers(1)));
    x.AB.KCa.add(strjoin(controllers(1)));
    x.AB.Kd.add(strjoin(controllers(1)));
    x.AB.NaV.add(strjoin(controllers(1)));
elseif controltype == 2
    x.AB.ACurrent.add(strjoin(controllers(2)));
    x.AB.CaS.add(strjoin(controllers(2)));
    x.AB.CaT.add(strjoin(controllers(2)));
    x.AB.HCurrent.add(strjoin(controllers(2)));
    x.AB.KCa.add(strjoin(controllers(2)));
    x.AB.Kd.add(strjoin(controllers(2)));
    x.AB.NaV.add(strjoin(controllers(2)));
else
    x.AB.ACurrent.add(strjoin(controllers(randi(2))));
    x.AB.CaS.add(strjoin(controllers(randi(2))));
    x.AB.CaT.add(strjoin(controllers(randi(2))));
    x.AB.HCurrent.add(strjoin(controllers(randi(2))));
    x.AB.KCa.add(strjoin(controllers(randi(2))));
    x.AB.Kd.add(strjoin(controllers(randi(2))));
    x.AB.NaV.add(strjoin(controllers(randi(2))));
end

% set controller variables
x.set('*tau_m',5e6./x.get('*gbar'));
x.set('*tau_filter',3e3);
x.set('*i_Ca_target',-119.8832); % see calcium_correlation.m

% add leak current to prinz model
x.AB.add('Leak','gbar',1e-4);

% set number of simulations to perform
simulationcount = 20;
conductances = zeros(simulationcount,8);

% run the simulations
x.t_end = 1e6;
for i = 1:simulationcount
    g0 = (2e-3/x.AB.A)*rand(8,1);
    x.set('*gbar',g0);
    g0 = (2e-3/x.AB.A)*rand(7,1);
    x.set('*Controller.m',g0);
    x.AB.Leak.gbar = (1e-4/x.AB.A)+((1.99e-2/x.AB.A)*rand());
    [~,~,M,~] = x.integrate;

    conductances(i,:) = x.get('*gbar');
    corelib.textbar(i,simulationcount);
end

% close all open plots
close all

% get a list of channel names and count them
channels = x.AB.find('conductance');
N = length(channels);

% plot channel correlations
figure(); hold on
ax = figlib.gridAxes(N);
c = lines;
idx = 1;
ph = gobjects(N-1,N);
for i = 1:N-1
    for j = i+1:N
        idx = idx + 1;
        ph(i,j) = scatter(ax(i,j),conductances(:,i),conductances(:,j),'MarkerFaceColor',c(idx,:),'MarkerEdgeColor',c(idx,:),'MarkerFaceAlpha',.2);
        if j < N
            set(ax(i,j),'XColor','w')
        end
        if i > 1
            set(ax(i,j),'YColor','w')
        end
        if j == N
            xlabel(ax(i,j),channels{i})
        end
        if i == 1
            ylabel(ax(i,j),channels{j})
        end
    end
end

figlib.pretty('PlotLineWidth',1)

drawnow

% plot convergence graph
figure('outerposition',[300 300 900 600],'PaperUnits','points','PaperSize',[1200 600]); hold on
subplot(2,1,1); hold on

time = x.dt*(1:length(M))*1e-3;
plot(time,M(:,2:2:end));
set(gca,'XScale','log','YScale','log','YTick',[1e-2 1e0 1e2 1e4])
xlabel('Time (s)')
ylabel('g (uS/mm^2)')

% plot trace
x.t_end = 5e3;
[V,~,~,~] = x.integrate;
subplot(2,1,2); hold on

time = x.dt*(1:length(V))*1e-3;
plot(time,V,'k')
set(gca,'YLim',[-80 50])
ylabel('V_m (mV)')
xlabel('Time (s)')

figlib.pretty('PlotLineWidth',1.5,'LineWidth',1.5)

drawnow
end