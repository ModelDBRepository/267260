clc
clear all
close all
%% Choose sample, current and neuron

names = ["1-female_footshock_first";"2-male_footshock_first";"3-female_group_first"; ...
    "4-female_single_first";"5-male_no-stress_first";"6-female_control_first";...
    "7-female_CORT_first";"8-male_control_first";"9-male_CORT_first"];

sample = 1;                               % choose sample from the names above: 1 to 9

load(names(sample))
parameters = firstbestparameters;         % load all fit neurons parameters
number_of_neurons = size(parameters,1);   % number of fit neuron in the file

neuron = 1;                               % select the neuron, from 1 to number_of_neurons
I = 30;                                   % selec the applied current
dt = 0.05;                                % Euler integration step

%% Extract experimental recording time and current for each sample

[T,in,time_sim]=exp_setup(sample,dt,I);

%% Run the Adex Model

[storeV,storeW,storeVt,storeVp,storeVr] = AdEx_integrator_Euler(parameters(neuron,:),T,dt,in);

%% plot model variables and current time evolution

fig = figure(2);   %plot training traces
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.9, 0.9])

subplot(3,2,1)
yyaxis left
plot(time_sim,storeV,'b','LineWidth',1.5), hold on
xlabel('$t$ (ms)','Interpreter','Latex','FontSize',18)
ylabel('$V$ (mV)','Interpreter','Latex','FontSize',18)
ylim([-150,50]);
yyaxis right
title(['Neuron ',num2str(neuron),'Sample ', names(sample)],'FontSize',15)
plot(time_sim,in,'r','LineWidth',1.5), hold on
xlim([0,T]);
ylim([-50,500]);
ylabel('$I$ (pA)','Interpreter','Latex','FontSize',18)


subplot(3,2,2)
plot(time_sim,storeW,'k','LineWidth',1.5), hold on
xlabel('$t$ (ms)','Interpreter','Latex','FontSize',18)
ylabel('$w$','Interpreter','Latex','FontSize',18)

subplot(3,2,3)
plot(time_sim,storeVt,'k','LineWidth',1.5), hold on
xlabel('$t$ (ms)','Interpreter','Latex','FontSize',18)
ylabel('$V_T$ (mV)','Interpreter','Latex','FontSize',18)

subplot(3,2,4)
plot(time_sim,storeVp,'k','LineWidth',1.5), hold on
xlabel('$t$ (ms)','Interpreter','Latex','FontSize',18)
ylabel('$V_p$ (mV)','Interpreter','Latex','FontSize',18)

subplot(3,2,5)
plot(time_sim,storeVr,'k','LineWidth',1.5), hold on
xlabel('$t$ (ms)','Interpreter','Latex','FontSize',18)
ylabel('$V_r$ (mV)','Interpreter','Latex','FontSize',18)


%    namefig = sprintf('V_w_Vt_Vp_Vr.svg');
%    print('-painters','-dsvg','-r300' ,'-loose',namefig)   %save a svg figure into the results folder

