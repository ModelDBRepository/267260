function [Vstore,current,timespike]=AdEx_integrator(par,N,nt,dt,transient,Irheo,freq,Amp,V0,w0,inhibition) 

I = zeros(N,1);                         % instantaneous mEPSC,
current= zeros(nt-3*transient/dt,1);
jumpstart = zeros(N,1);                 % save last EPSC spike time
Amplitude = Amp;
Iraux=zeros(N,1);
Vstore = zeros(nt-3*transient/dt,1);
timespike = zeros(nt*10,2);             % save neuron index and spike time
transient=transient/dt;

C = par(:,1);                           % capacitance
gl = par(:,2);                          % leak conductance
El = par(:,3);                          % leak reversal potential
Vt0 = par(:,4);                         % resting spike threshold
delta_t = par(:,5);                     % slope factor
tau_w = par(:,6);                       % adaptation time constant
a = par(:,7);                           % subthreshold adaptation
b = par(:,8);                           % spike-triggered adaptation
Vr0 = par(:,9);                         % linear coefficient of the resting reset voltage : Vr0 = Vra + Vrb*I
tau_vt = par(:,10);                     % time constant for Vt    
q = par(:,11);                          % reset for Vt
Vp0 = par(:,12);                        % linear coefficient of the peak voltage : Vp0 = Vpa + Vpb*I
tau_p=par(:,13);                        % max peak time constant 
p=par(:,14);                            % max peak adaptation
tau_r=par(:,15);                        % reset voltage time constant
r = par(:,16);                          % reset voltage adaptation

Vt0 = Vt0+Vr0;                          % resting spike threshold

%%%%%%%%%%%%%%%%%%initial conditions

V = V0;                                 % initial voltage
w = w0;                                 % initial adaptation variable value
Vt = Vt0;                               % initial threshold
Vp =  Vp0;                              % initial peak voltage for reset
Vr = Vr0;                               % initial reset voltage

%% Integration start

ns=0;
chave=0;
cont=0;    
cont3=0;                                % auxiliar variable to print time at the screen

for i = 1:nt   
%---------------- generate mEPSC  

    if (i>transient && i<=3*transient)
        I=Irheo;
    elseif i>3*transient

        cont=cont+1;

        qq = rand(N,1);                            % random number for each neuron between 0-1, to determine which neuron will receive a current pulse
        aa = rand(N,1);                            % random number for each neuron between 0-1, to define if the pulse will be excitatory or inhibitory
        prob = freq*(dt/1000);                     % fixed probability of receiveing a current pulse
        jump = find(qq<prob);                      % neurons that will receive a pulse
        if size(jump,1)>0
            jumpstart(jump)=dt*i;                   % time start of the jump, updated every time a neuron receives a new pulse
            Amplitude(jump) = Amp(jump).*((aa(jump)>=inhibition).*1)-Amp(jump).*((aa(jump)<inhibition).*1); %define the amplitude being excitatory or inhibitory
        end

        Itd=2.0;                                   % decay time of the current pulse
        Itr=0.5;                                   % rise time ofthe current pulse
        t_norm=((Itd*Itr)/(Itd-Itr))*log(Itd/Itr); % normalization factor
        Inorm_factor = (exp(-t_norm./Itd)-exp(-t_norm./Itr));  %normalization factor

        Iraux = Iraux.*exp(-dt./Itd)+(I-Irheo).*((i*dt-jumpstart)==0);  % auxiliar to store the previous current, before the arrivel of a new pulse
        I=Irheo+Iraux+(Amplitude./Inorm_factor).*(exp(-(i*dt-jumpstart)./Itd)-exp(-(i*dt-jumpstart)./Itr));  % input current

        Vstore(cont)=V(N);                         % store voltage of neuron N
        current(cont)=I(N);                        % store the input current of neuron N
    end

    %%----------------------------
    %%%%%%%%%%%%%%%%%%%%%% Euler Integrator - Adex + dVt +dVp +dVr
    V = V + dt*(-gl.*(V-El)+gl.*delta_t.*exp((V-Vt)./delta_t)-w+I)./C;
    w = w + dt*(a.*(V-El)-w)./tau_w;

    Vt= Vt+ dt*((Vt0-Vt)./tau_vt);
    Vp= Vp+ dt*((Vp0-Vp)./tau_p);
    Vr= Vr+ dt*((Vr0-Vr)./tau_r);
    %%%%%%%%%%%%%%%%%%%%%%

    if i>3*transient
        spikes = find(V>=Vp);                                              %find neurons that spike in that time instant
        if ~isempty(spikes)==1
            timespike(ns+1:ns+length(spikes),:) = [spikes,spikes*0+dt*i];  %save neuron index and spike time, for raster plot
            ns = ns + length(spikes);                                      %total number of spikes
        end
    end


    chave=chave+1;
    if chave == 100000                  % print time at the screen
        chave=0;
        cont3 = cont3+1000;
        fprintf('%d ms\n',cont3)
    end

    % update variables when voltage reaches the peak limit
    w = w + b.*(V>=Vp);
    Vt = Vt + q.*(V>=Vp);
    V2= V;
    V = V + (Vr-V).*(V>=Vp);
    Vp = Vp + p.*(V2>=Vp).*(Vp>Vt);
    Vp = Vp - p.*(V2>=Vp).*(Vp<Vt);      % avoid that Vp becames smaller than Vt
    Vr = Vr + r.*(V2>=Vp);
end
%% 

timespike = timespike(1:ns,:);        % rescale the size of the vector to free memory space
end