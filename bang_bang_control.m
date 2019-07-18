
close all

% reproduces Fig 3c in Tim O'Leary's paper

x = xolotl.examples.BurstingNeuron('prefix','liu');
    
x.AB.Ca_target = 7;
    
x.AB.NaV.add('plausible-homeostasis/BangBangController','tau_m',666);
x.AB.CaT.add('plausible-homeostasis/BangBangController','tau_m',55555);
x.AB.CaS.add('plausible-homeostasis/BangBangController','tau_m',45454);
x.AB.ACurrent.add('plausible-homeostasis/BangBangController','tau_m',5000);
x.AB.KCa.add('plausible-homeostasis/BangBangController','tau_m',1250);
x.AB.Kd.add('plausible-homeostasis/BangBangController','tau_m',2000);
x.AB.HCurrent.add('plausible-homeostasis/BangBangController','tau_m',125000);

x.t_end = 5e5;
x.sim_dt = .1;
x.dt = 100;

simulationcount = 100;
conductances = zeros(simulationcount,8);

for i = 1:simulationcount
    g0 = (0.2/x.AB.A)*rand(8,1);
    x.set('*gbar',g0);
    x.set('*Controller.m',0);
    x.AB.Leak.gbar = (0.001/x.AB.A)+((0.199/x.AB.A)*rand());

    [~,~,C] = x.integrate;

    conductances(i,:) = x.get('*gbar');
    corelib.textbar(i,simulationcount);
end

conductances(:,7) = [];

figure(); hold on
plotmatrix(conductances)

drawnow

figure('outerposition',[300 300 900 600],'PaperUnits','points','PaperSize',[1200 600]); hold on
subplot(2,1,1); hold on

time = x.dt*(1:length(C))*1e-3;
plot(time,C(:,2:2:end));
set(gca,'XScale','log','YScale','log','YTick',[1e-2 1e0 1e2 1e4])
xlabel('Time (s)')
ylabel('g (uS/mm^2)')

subplot(2,1,2); hold on
x.dt = .1;
x.t_end = 1e3;
V = x.integrate;
time = x.dt*(1:length(V))*1e-3;
plot(time,V,'k')
set(gca,'YLim',[-80 50])
ylabel('V_m (mV)')
xlabel('Time (s)')

drawnow

figlib.pretty('PlotLineWidth',1.5,'LineWidth',1.5)