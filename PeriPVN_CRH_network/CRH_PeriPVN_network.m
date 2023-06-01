clear all;close all;clc
%%
tic

load('Irheo_analytical_and_bifurc.mat')     % load file with rheobases(first column, in pA) and bifurcation type (second columnn 1 = Saddle, 0 = Hopf)
                                            % rows 1-151 Protocol 1 neurons, rows 152-248 Protocol 2 neurons
load('parameters.mat')                      % load fit parameters of all neurons
                                            % column order (C, gl, El, vt0, delta, tau_w, a, b, Vr0, tau_vt, q, Vp0, tau_p, p, tau_r, r, fit error)
                                            % rows 1-151 Protocol 1 neurons, rows 152-248 Protocol 2 neurons


results = sprintf('results/');              % create folder so save results
mkdir(results);                             % create the directory

%% Network Parameters and Initial conditions

N = 200;                                    % total number of neurons, each network = N/2
model_aux = [ones(N/2,1);zeros(N/2,1)];     % auxiliar to specify neuron model, equals 1 for the modified AdEx and 0 for LIF
cp = 0.1;                                   % coupling probability
cw = 2;                                     % coupling weight
Hopf_prop = 0.4;                            % proportion probability of Hopf neuron in the CRH network

Iext = zeros(N,1);                          % external current vector
Iext(N/2+1:N,1) = 20;                       % baseline current for the PeriPVN interneurons
A = 100;                                    % amplitude of the postsynaptic current in pA

%% Coupling Matrix
Wij = (rand(N,N)<=cp).*(1/(cp*sqrt(N)))*cw; % random NxN matrix
Wij(1:N/2,1:N/2) = 0;                       % uncounpling the index inside CRH network
Wij(:,N/2+1:end)=-Wij(:,N/2+1:end);         % turns all interneurons connections negative

%% Integration parameters

dt = 0.01;                                  % integraton step
T = 10000;                                  % simulation time in ms
nt = round(T/dt);                           % total integration steps

step=1;                                     % time step to store V
time_sim = (1:nt)*dt;                       % vector with all sim steps
%% Selection of parameters for each neuron

par = zeros(N/2,16);                        % defining the CRH parameters matrix
bif = zeros(N/2,1);                         % CRH bifurcation type vector

for i = 1:(1-Hopf_prop)*N/2                 % selection of parameters for the Saddle neurons in the CRH network
    aux =round(150*rand())+1;                   % auxiliar variable in the range 1-151, to select only Protocol 1 neurons
    while Irheo_bif(aux,2)==0               % random select a neuron until it is a SADDLE, Irheo_bif(aux,2)=1
        aux =round(150*rand())+1;
    end
    par(i,:)=parameters(aux,1:16);          % save parameters
    bif(i)= Irheo_bif(aux,2);               % 1 = Saddle node bif, 0 = Hopf bif.
    Iext(i)=Irheo_bif(aux,1)-0.1;           % save respective rheobase
end

for i = (1-Hopf_prop)*N/2+1:N/2             % selection of parameters for the Hopf neurons in the CRH network
    aux =round(150*rand())+1;               % % auxiliar variable in the range 1-151, to select only Protocol 1 neurons
    while Irheo_bif(aux,2)==1               % random select a neuron until it is a Hopf, Irheo_bif(aux,2)=0
        aux =round(150*rand())+1;
    end
    par(i,:)=parameters(aux,1:16);          % save parameters
    bif(i)=  Irheo_bif(aux,2);              % 1 = Saddle node bif, 0 = Hopf bif.
    Iext(i)=Irheo_bif(aux,1)-0.1;           % save respective rheobase
end

[hopfs,~]=find(bif==0);                     % find indexes of all hopfs neurons in bif
[saddles,~]=find(bif==1);                   % find indexes of all saddles neurons in bif

%% separating saddle from hopf

par_saddle=par(saddles,:);                  % select only the saddles in par
par_hopf = par(hopfs,:);                    % select only the hopfs in par

%% Peri_PVN neurons  - LIF model
% C  gl   El            Vr      Vp
par_PeriPVN = [20 0.3 -70 1 1 1 1 1 -60 1 1 -20 1 1 1 1].*ones(N/2,16);  % PeriPVN neuron LIF parameters,


%%
par_new = [par_saddle;par_hopf;par_PeriPVN];   % all neurons parameters sorted by type Saddle-Hopf-Interneuron

%% file names

namefig = sprintf('%sRaster_CW%.1f_PSP_Amplitude%.1f_Hopf_Probability%.1f.svg',results,cw,A, Hopf_prop);
namefile= sprintf('%sRaster_CW%.1f_PSP_Amplitude%.1f_Hopf_Probability%.1f.mat',results,cw,A, Hopf_prop);
text1 = sprintf('Raster, cw = %.1f pA, PSP Amplitude = %.1f pA, Hopf Probability = %.1f',cw,A, Hopf_prop);
%% Simulation

[timespike]=AdEx_integrator_CRH_PeriPVN(par_new,N,nt,dt,Iext,A,model_aux,Wij);

%% Saving raster file

raster=[timespike(:,2),timespike(:,1)];          % time in ms, neuron index
save(namefile, 'raster');

%% raster plot and PSTH of the CRH network

[saddle_idx,~]=find(raster(:,2)<=length(saddles));
[hopf_idx,~]= find(raster(:,2)>length(saddles) & raster(:,2)<=N/2);
[Peri_PVN_idx,~]=find(raster(:,2)>N/2);
[CRH_idx,~]=find(raster(:,2)<=N/2);

raster_saddle = raster(saddle_idx,:);
raster_hopf= raster(hopf_idx,:);
raster_Peri_PVN = raster(Peri_PVN_idx,:);

fig1 = figure(1);
subplot(11,1,[1 2 3 4 5])
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.05, 1.0, 0.9])
plot(raster_saddle(1:1:end,1),raster_saddle(1:1:end,2),'.','Markersize',10,'MarkerEdgeColor',[0.8500 0.3250 0.0980],'MarkerFaceColor',[0.8500 0.3250 0.0980]),hold on
plot(raster_hopf(1:1:end,1),raster_hopf(1:1:end,2),'.','Markersize',10,'MarkerEdgeColor',[0 0.4470 0.7410],'MarkerFaceColor',[0.9290 0.6940 0.1250]),hold on
plot(raster_Peri_PVN(1:5:end,1),raster_Peri_PVN(1:5:end,2),'.','Markersize',10,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0]),hold on
set(gca,'TickLabelInterpreter','latex','FontSize',20)
legend('Saddle','Hopf','PeriPVN','Interpreter','Latex','FontSize',15)
legend('boxoff')
ylabel('Neuron','Interpreter','Latex','FontSize',20)
xlabel('time (ms)','Interpreter','Latex','FontSize',20)
ylim([0 N+0.4*N])
xlim([0 T])
title(text1,'Interpreter','Latex','FontSize',20)

bin=50;
subplot(11,1,[8 9 10 11])
[h,x] = hist(raster(CRH_idx,1),min(raster(CRH_idx,1)):bin:max(raster(CRH_idx,1)));
h = (1000*h)/((N/2)*bin);
plot(x,h,'-','LineWidth',2,'Color',[0 0.4470 0.7410]), hold on
set(gca,'TickLabelInterpreter','latex','FontSize',20)
ylabel('PSTH (Hz)','Interpreter','Latex','FontSize',20) %(x$10^{-4}$)
xlabel('time (ms)','Interpreter','Latex','FontSize',20)
ylim([-0.2 2.5])


print('-painters','-dsvg','-r300' ,'-loose',namefig);   %save a png figure into the results folde
% clf(fig1)


%%
toc
