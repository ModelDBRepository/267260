function [storeV,storeW,storeVt,storeVp,storeVr]=AdEx_integrator_Euler(par,T,dt,in)


nt = round(T/dt);   %number of iteractions

C = par(1);           %capacitance
gl = par(2);          %leak conductance
El = par(3);          %leak reversal potential
Vt0 = par(4);         %resting spike threshold
delta_t = par(5);     %slope factor
tau_w = par(6);       %adaptation time constant
a = par(7);           %subthreshold adaptation
b = par(8);           %spike-triggered adaptation
Vr0 = par(9);         %resting reset voltage
tau_vt = par(10);     %time constant for Vt
q = par(11);          %reset for Vt
Vp0 = par(12);        %resting max peak voltage
tau_p=par(13);        %max peak time constant
p=par(14);            %max peak adaptation
tau_r=par(15);        %reset voltage time constant
r = par(16);          %reset voltage adaptation

Vt0 = Vt0+Vr0;        %adjusting resting spike threshold to be bigger than the reset

%initial conditions
V = -50;
w = 0;
Vt = Vt0;
Vp = Vp0;
Vr = Vr0;

%store variables
storeV = zeros(1,nt);
storeW = zeros(1,nt);
storeVt = zeros(1,nt);
storeVp = zeros(1,nt);
storeVr = zeros(1,nt);
%% integrate and store the AdEx + dVt
for i = 1:nt
    
    I = in(:,i);
    
    %----------------- Euler Integrator - Adex + dV_T
    
    V = V + dt*(-gl.*(V-El)+gl.*delta_t.*exp((V-Vt)./delta_t)-w+I)./C;
    w = w + dt*(a.*(V-El)-w)./tau_w;
    Vt= Vt+ dt*((Vt0-Vt)./tau_vt);
    
    Vp= Vp+ dt*((Vp0-Vp)./tau_p);
    Vr= Vr+ dt*((Vr0-Vr)./tau_r);
    
    %--------------------------------------
    
    Vs = V + (Vp-V).*(V>=Vp);  %adjust the max voltage to the peak voltage
    
    %-------------------------------------
    %store variables
    storeV(:,i) = Vs;          %save the voltage v time evolution
    storeW(:,i)= w;            %save w time evolution
    storeVt(:,i)= Vt;          %save VT time evolution
    storeVp(:,i)= Vp;          %save Vp time evolution
    storeVr(:,i)= Vr;          %save Vr time evolution
    
    % updates after spike
    w = w + b.*(V>=Vp);
    Vt = Vt + q.*(V>=Vp);
    V2= V;
    V = V + (Vr-V).*(V>=Vp);
    Vp = Vp + p.*(V2>=Vp);
    Vr = Vr + r.*(V2>=Vp);
end

end