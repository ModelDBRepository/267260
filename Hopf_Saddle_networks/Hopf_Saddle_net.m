clear all;close all;clc
%%
tic

load('Irheo_analytical_and_bifurc.mat')              % load file with rheobases(first column, in pA) and bifurcation type (second columnn 1 = Saddle, 0 = Hopf)
                                                     % rows 1-151 Protocol 1 neurons, rows 152-248 Protocol 2 neurons 
load('parameters.mat')                               % load fit parameters of all neurons
                                                     % column order (C, gl, El, vt0, delta, tau_w, a, b, Vr0, tau_vt, q, Vp0, tau_p, p, tau_r, r, fit error)                                                        
                                                     % rows 1-151 Protocol 1 neurons, rows 152-248 Protocol 2 neurons 

load('stable_fixed_points_bif_freq-analytical.mat')  % load neurons fixed points V and w for the initial conditions

%% create directory to save files and plots

results = sprintf('results/');                       % create folder so save results
mkdir(results);                                      % create the directory

%% 
sim = 0;                                             % auxiliar variable
inhi_steps = 5;                                      % amount of point in the interval 0-1 of the inhibition proportion
mean_activity = zeros(1+inhi_steps,1);               % mean activity for each inhibition proportion
for inhibition = 0:1/inhi_steps:1                    % loop over the inhibition proportion 
    %% Network Parameters
    sim = sim +1;
    N = 50;                                          % number of neuron in the network

    par = zeros(N,16);                               % store neurons parameters
    bif = zeros(N,1);                                % store bifurcation type Saddle = 1, Hopf = 0
    Irheo = zeros(N,1);                              % store rheobase current for each neuron
    V0=zeros(N,1);                                   % initial condition of V
    w0=zeros(N,1);                                   % initial condition of w

    neuron_type = 0;                                 % Hopf = 1,   Saddle = 0 

    %% selecting neurons and their respective parameters

    for i = 1:N
        aux =round(150*rand())+1;                            % auxiliar variable in the range 1-151, selects neurons from Protocol 1 only
            while stable_v_w_bif_freq(aux,3)==neuron_type    
                aux =round(150*rand())+1; 
            end   
        par(i,:)=parameters(aux,1:16);                       % select parameter from PSO fit neurons
        bif(i)=stable_v_w_bif_freq(aux,3);                   % 1 = Saddle node bif, 0 = Hopf bif.
        Irheo(i)=Irheo_bif(aux,1)-0.2;                       % select baseline current of each neuron, below rheobase 
        V0(i)=stable_v_w_bif_freq(aux,1);                    % v of the fixed point 
        w0(i)=stable_v_w_bif_freq(aux,2);                    % w of the fixed point
    end

    %% Integration parameters

    dt = 0.01;                                           % integraton step
    transient = 10000;                                   % transient to start appling current
    T = 20000;                                           % simulation time in ms
    nt = round((T+3*transient)/dt);                      % total integration steps
    time_sim = (1:nt-3*transient/dt)*dt;                 % vector with all sim steps

    %% I parameters

    freq=10;                                             % mean frequency of pulses in Hz
    Amp=30.*ones(N,1);                                   % Amplitude of current pulses in pA

    %% Simulation

    [V,current,timespike]=AdEx_integrator(par,N,nt,dt,transient,Irheo,freq,Amp,V0,w0,inhibition);

    %% Saving files

    raster=[timespike(:,2),timespike(:,1)];
    namefile= sprintf('%sraster_inhibition%.2f.mat',results,inhibition);
    save(namefile, 'raster'); 

    %% mean activity

    number_spikes = find(timespike(:,2)>(nt*dt-10000));    % calculate the number of spikes in the last 10s of simulationn
    mean_activity(sim) = length(number_spikes)/(10*N);     % calculates the mean activity

    %% raster plot and PSTH plot

    fig1 = figure(1);
    subplot(11,1,[1 2 3 4 5])
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.05, 1.0, 0.9])

    text1 = sprintf('I:  Amplitude = %.0d pA, freq = %.0d Hz, %.2f exc, %.2f inhi ',Amp(1),freq,1-inhibition,inhibition);
    plot(raster(1:1:end,1)-3*transient,raster(1:1:end,2),'.','Markersize',10,'MarkerEdgeColor',[0.8500 0.3250 0.0980],'MarkerFaceColor',[0.8500 0.3250 0.0980]),hold on
    set(gca,'TickLabelInterpreter','latex','FontSize',20)
    legend('Saddle','Hopf','Interpreter','Latex','FontSize',15)
    legend('boxoff')
    ylabel('Neuron','Interpreter','Latex','FontSize',20)
    xlabel('time (ms)','Interpreter','Latex','FontSize',20)
    ylim([0 N+50])
    title(text1,'Interpreter','Latex','FontSize',20)

    bin=100;
    subplot(11,1,[7 8 9 10 11])
    [h,x] = hist(raster(:,1),min(raster(:,1)):bin:max(raster(:,1)));
    h = (1000*h)/(N*bin);
    plot(x-3*transient,h,'-','LineWidth',2,'Color',[0 0.4470 0.7410]), hold on 
    set(gca,'TickLabelInterpreter','latex','FontSize',20)
    ylabel('PSTH (Hz)','Interpreter','Latex','FontSize',20) %(x$10^{-4}$)
    xlabel('time (ms)','Interpreter','Latex','FontSize',20)

    namefig = sprintf('%sraster_PSTH_inhiprob_%.2f.svg',results,inhibition);
    print('-painters','-dsvg','-r300' ,'-loose',namefig);   %save a png figure into the results folde
    clf(fig1)

    %% Voltage trace and inut current of a random neuron

    fig2 = figure(2);
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.05, 0.9, 0.9])
    subplot(2,1,1)
    plot(time_sim,V,'b-','Linewidth',1.5)
    xlabel('$time$ (ms)','Interpreter','Latex','FontSize',20)
    ylabel('$V$ (mV)','Interpreter','Latex','FontSize',20)
    subplot(2,1,2)
    plot(time_sim,current,'r-','Linewidth',1.5)
    xlabel('$time$ (ms)','Interpreter','Latex','FontSize',20)
    ylabel('$I$ (pA)','Interpreter','Latex','FontSize',20)

    namefig3 = sprintf('%sNeuron_trace_inputcurrent_inhiprob_%.1d.svg',results,inhibition);
    print('-painters','-dsvg','-r300' ,'-loose',namefig3);   %save a png figure into the results folde
    clf(fig2)

end
%% Mean activity plot

proportion = (0:0.2:1);
figure(3)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.05, 0.9, 0.9])
plot(proportion,mean_activity,'r-o','Linewidth',4)
xlabel('Inhibition Proportion','Interpreter','Latex','FontSize',20)
ylabel('Mean Activity (Hz)','Interpreter','Latex','FontSize',20)
ylim([0 2.8])

namefig2= sprintf('%smean_activity.svg',results);
print('-painters','-dsvg','-r300' ,'-loose',namefig2);   %save a png figure into the results folde
inhi_meanactivity= [proportion',mean_activity];
namefile2= sprintf('%sinhibition_mean-activity.mat',results);
save(namefile2, 'inhi_meanactivity'); 


%%
toc
