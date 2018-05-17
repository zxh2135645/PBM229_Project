% First version of IR simulation
clear all;
close all;
clc;

df = 0;		% Hz off-resonance.
T1 = 600;	% ms.
T2 = 100;	% ms.
TR = 2000;     % ms.
TI = 750;
TE = 10;		% ms.
dt = 0.5;   % ms.
flip90 = pi/2;	% radians.
flip180 = pi; % radians


Mz0 = 1;
M = [0;0;Mz0];
[Adt,Bdt] = freeprecess(dt,T1,T2,df);
step = 1;
NTR = 1; % number of TRs to sim
NstepTR = TR / dt; % sim steps per TR
NstepTE = TE / dt; 
NstepTI = TI / dt;
NstepTot= NTR * NstepTR; % total sim steps
Msim = zeros(3, NstepTot); %keep track of temporal evolution

% Sim TR by TR
% assuming instantaneous RF excitation ?and grad spoiling?
for nn = 1: NTR
    Mnew = yrot(flip180) * M; % 180 degrees
    Mnew = Adt * Mnew + Bdt;
    Msim(:,step) = Mnew; step = step + 1;
    
    % Inversion
    for ss = 2:(NstepTI-1)
        Mnew = Adt*Msim(:, step-1) + Bdt;
        Msim(:,step) = Mnew; step = step + 1;
    end
    
    Mnew = yrot(flip90) * Msim(:, step-1); % 90 degrees
    
    % first free precession step in this TR 
    % included in the same time interval as excitation
    Mnew = Adt * Mnew + Bdt;
    Msim(:,step) = Mnew; step = step + 1;
    
    % subsequent steps include free precession
    for ss = (NstepTI+1):(NstepTI+NstepTE/2-1)
        Mnew = Adt*Msim(:, step-1) + Bdt;
        % final step includes gradient spoiling?
        
        Msim(:,step) = Mnew; step = step + 1;
    end
    
    Mnew = yrot(flip180) * Msim(:, step-1); % 180 degrees
    Mnew = Adt * Mnew + Bdt;
    Msim(:,step) = Mnew; step = step + 1;
    
    for ss = (NstepTI+NstepTE/2+1): NstepTR
        Mnew = Adt*Msim(:, step-1) + Bdt;
        % final step includes gradient spoiling?
        
        Msim(:,step) = Mnew; step = step + 1;
    end
end
Mxy = Msim(1,:) + 1i*Msim(2,:);

time = (1:NstepTot) * dt;
figure(); plot( time,Msim(1,:), time,Msim(2,:), time, Mxy, time,Msim(3,:), ...
    'LineWidth', 1);
legend('Mx', 'My', 'Mxy', 'Mz'); ylim([-1 1]); grid on;
title('Saturation Recovery');
xlabel('ms'); ylabel('Normalized Magnetization');