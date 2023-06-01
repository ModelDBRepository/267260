function [T,in,time_sim]=exp_setup(sample,dt,I)


%% female footshock sample
if sample==1
    
    tin1 = 550;                             % first current injection start
    tout1 = tin1 + 500;                     % first current injection stop
    tin2 = 1050;                            % second current injection start
    tout2 = tin2 + 500;                     % second current injection stop
    
    T = 3000;                               % max time
    nt = round(T/dt);                       % total integration time for sim
    time_sim = (1:nt)*dt;                   % vector with all time steps
    
    in = (time_sim>tin1).*(time_sim<tout1)*(-30)+(time_sim>tin2).*(time_sim<tout2)*I; % current vector for all simulation steps
end

%% male footshock sample
if sample==2
    
    tin1 = 250;                             % first current injection start
    tout1 = tin1 + 400;                     % first current injection stop
    tin2 = 650;                             % second current injection start
    tout2 = tin2 + 500;                     % second current injection stop
    
    T = 3000;                               % max time
    nt = round(T/dt);                       % total integration time for sim
    time_sim = (1:nt)*dt;                   % vector with all time steps
    
    in = (time_sim>tin1).*(time_sim<tout1)*(-30)+(time_sim>tin2).*(time_sim<tout2)*I; % current vector for all simulation steps
end

%% female group and female single samples
if (sample==3 || sample==4)
    
    tin1 = 150;                             % first current injection start
    tout1 = tin1 + 400;                     % first current injection stop
    tin2 = 550;                             % Second current injection start
    tout2 = tin2 + 500;                     % Second current injection stop
    
    T = 3000;                               % max time
    nt = round(T/dt);                       % total integration time for sim
    time_sim = (1:nt)*dt;                   % vector with all time steps
    
    in = (time_sim>tin1).*(time_sim<tout1)*(-30)+(time_sim>tin2).*(time_sim<tout2)*I; % current vector for all simulation steps
end
%% Male non stress, Female and Male CORT, Female and Male Control samples
if sample>=5
    
    tin = 115;                             % current injection start
    tout = tin + 500;                      % current injection stop
    
    T = 1000;                              % max time
    nt = round(T/dt);                      % total integration time for sim
    time_sim = (1:nt)*dt;                  % vector with all time steps
    
    in = (time_sim>tin).*(time_sim<tout)*I;% current vector for all simulation steps
end


end
