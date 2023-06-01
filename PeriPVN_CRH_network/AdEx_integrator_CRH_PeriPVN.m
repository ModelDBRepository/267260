function [timespike]=AdEx_integrator_CRH_PeriPVN(par,N,nt,dt,Iext,A,model_aux,Wij)


PSP = zeros(N,1);             % PSP current os each neuron
raux=zeros(N,1);              % auxiliar variable for the PSC currents
lastspike=zeros(N,1);         % variable store the last EPSC spike
timespike = zeros(nt*10,2);   % variable to store neuron index and spike time

C = par(:,1);                 % capacitance
gl = par(:,2);                % leak conductance
El = par(:,3);                % leak reversal potential
Vt0 = par(:,4);               % resting spike threshold
delta_t = par(:,5);           % slope factor
tau_w = par(:,6);             % adaptation time constant
a = par(:,7);                 % subthreshold adaptation
b = par(:,8);                 % spike-triggered adaptation
Vr0 = par(:,9);               % linear coefficient of the resting reset voltage : Vr0 = Vra + Vrb*I
tau_vt = par(:,10);           % time constant for Vt
q = par(:,11);                % reset for Vt
Vp0 = par(:,12);              % linear coefficient of the peak voltage : Vp0 = Vpa + Vpb*I
tau_p=par(:,13);              % max peak time constant
p=par(:,14);                  % max peak adaptation
tau_r=par(:,15);              % reset voltage time constant
r = par(:,16);                % reset voltage adaptation


Vt0 = Vt0+Vr0.*model_aux;     % resting spike threshold

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%initial conditions

V = -60+10*rand(N,1);
w = zeros(N,1);
Vt = Vt0;
Vp =  Vp0;
Vr = Vr0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Integration start

ns=0;
chave=0;                    % auxiliar variable to print time at the screen
cont=0;                     % auxiliar variable to print time at the screen

for i = 1:nt                % time loop
    
    all_PSP = Wij*PSP;      % sum of all PSP inputs
    
    %%----------------------------
    %%%%%%%%%%%%%%%%%%%%%% Euler Integrator - Adex + dVt +dVp +dVr
    V = V + dt*(-gl.*(V-El)+model_aux.*(gl.*delta_t.*exp((V-Vt)./delta_t)-w)+Iext+all_PSP)./C;
    w = w + dt*(a.*(V-El)-w)./tau_w;
    
    Vt= Vt+ dt*model_aux.*((Vt0-Vt)./tau_vt);
    Vp= Vp+ dt*model_aux.*((Vp0-Vp)./tau_p);
    Vr= Vr+ dt*model_aux.*((Vr0-Vr)./tau_r);
    %%%%%%%%%%%%%%%%%%%%%%
    
    spikes = find(V>=Vp);                                              %find neurons that spike in that time instant
    if ~isempty(spikes)==1
        timespike(ns+1:ns+length(spikes),:) = [spikes,spikes*0+dt*i];  %save neuron index and spike time, for raster plot
        ns = ns + length(spikes);                                      %total number of spikes
        lastspike(spikes)=dt*i;                                        % save the last EPSC spike
    end
    
    %%%%%%PSCs of all neurons
    td=2.0;                                                            % decay time im ms
    tr=1.0;                                                            % rise time in ms
    t_norm2=((td*tr)/(td-tr))*log(td/tr);                              % normalization factor
    norm_factor = (exp(-t_norm2./td)-exp(-t_norm2./tr));               % normalization factor
    
    raux = raux.*exp(-dt./td)+(PSP-raux).*((i*dt-lastspike)==0);       % save the PSC value before the new spike
    PSP=raux+ (A./norm_factor).*(exp(-(i*dt-lastspike)./td)-exp(-(i*dt-lastspike)./tr)).*(lastspike>0); % PSC of each neuron
    
    PSP(isnan(PSP))=0;
    %%%%%%%%%%%
    
    
    chave=chave+1;
    if chave == 100000                                                % print time at the screen every 1000ms
        chave=0;
        cont = cont+1000;
        fprintf('%d ms\n',cont)
    end
    
    % update variables when voltage reaches the peak limit
    w = w + b.*(V>=Vp);
    Vt = Vt + model_aux.*q.*(V>=Vp);
    V2= V;
    V = V + (Vr-V).*(V>=Vp);
    Vp = Vp + model_aux.*p.*(V2>=Vp).*(Vp>Vt);
    Vp = Vp - model_aux.*p.*(V2>=Vp).*(Vp<Vt);                         % avoid that Vp becames smaller than Vt
    Vr = Vr + model_aux.*r.*(V2>=Vp);
    
    
end
%%

timespike = timespike(1:ns,:);   % rescale the size of the vector to free memory space
end