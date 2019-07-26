
close all

% reproduces Fig 3c in Tim O'Leary's paper

x = xolotl.examples.BurstingNeuron('prefix','prinz');
x.AB.add('Leak','gbar',1e-5);

x.t_end = 5e3;
[V,Ca,M,I] = x.integrate;
x.AB.Ca_target = x.AB.Ca_average;
    
x.AB.NaV.add('oleary/IntegralController','tau_m',5e3/x.AB.NaV.gbar);
x.AB.CaT.add('oleary/IntegralController','tau_m',5e3/x.AB.CaT.gbar);
x.AB.CaS.add('oleary/IntegralController','tau_m',5e3/x.AB.CaS.gbar);
x.AB.ACurrent.add('oleary/IntegralController','tau_m',5e3/x.AB.ACurrent.gbar);
x.AB.KCa.add('oleary/IntegralController','tau_m',5e3/x.AB.KCa.gbar);
x.AB.Kd.add('oleary/IntegralController','tau_m',5e3/x.AB.Kd.gbar);
x.AB.HCurrent.add('oleary/IntegralController','tau_m',5e3/x.AB.HCurrent.gbar);

x.t_end = 5e5;
x.sim_dt = .1;
x.dt = 100;

simulationcount = 10;
conductances = zeros(simulationcount,8);

for i = 1:simulationcount
    g0 = (2e-3/x.AB.A)*rand(8,1);
    x.set('*gbar',g0);
    x.set('*Controller.m',0);
    x.AB.Leak.gbar = (1e-5/x.AB.A)+((1.99e-3/x.AB.A)*rand());
    [V,Ca,M,I] = x.integrate;

    conductances(i,:) = x.get('*gbar');
    corelib.textbar(i,simulationcount);
end

conductances(:,7) = [];
figure(); hold on
plotmatrix(conductances)

drawnow

figure('outerposition',[300 300 900 600],'PaperUnits','points','PaperSize',[1200 600]); hold on
subplot(2,1,1); hold on

time = x.dt*(1:length(M))*1e-3;
plot(time,M(:,2:2:end));
set(gca,'XScale','log','YScale','log','YTick',[1e-2 1e0 1e2 1e4])
xlabel('Time (s)')
ylabel('g (uS/mm^2)')

subplot(2,1,2); hold on
x.dt = .1;
x.t_end = 5e3;
V = x.integrate;
time = x.dt*(1:length(V))*1e-3;
plot(time,V,'k')
set(gca,'YLim',[-80 50])
ylabel('V_m (mV)')
xlabel('Time (s)')

drawnow

figlib.pretty('PlotLineWidth',1.5,'LineWidth',1.5)