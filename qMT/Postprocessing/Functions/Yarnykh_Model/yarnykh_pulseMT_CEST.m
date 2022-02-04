function [Mz,Mzn,Mz0] = yarnykh_pulseMT_CEST(M0,R1,T2,TR,kba,tm,ts,thetaEX,deltaMT,B1e,lineshape)
%% This function is designed to find Magnetization constants from a set of
% data for a qMT image analysis. From Yarnykh, 2004
% Author: Alex Smith
% Date: 2/14/13
% Cleaning up Code - 2/14/14

% B1e - RF Power (uT)
% thetaEX - Excitation Angle, radians

% Make sure arrays are proper dim
if size(M0,1) == 1
    M0 = M0';
end

% CEST offsets
deltaCEST = -640:64:640;
B1e = [B1e,2];

% Set proton gyromagenetic ratio (rad/s/uT)
gamma = 42.58*2*pi;


% Calc RF saturation rate
% Liquid Pool
ga = absorptionLineShape(T2(1),[deltaCEST,deltaMT],'lorentzian');
Wa = pi*gamma^2*ga*B1e.^2;
% Wa = zeros(size(ga));

% Bound Pool
gb = absorptionLineShape(T2(2),deltaMT,lineshape);
Wb = pi*gamma^2*gb*B1e.^2;

%extend to CEST range
gb2 = absorptionLineShape(T2(2),[deltaCEST,deltaMT],lineshape);
Wb2 = pi*gamma^2*gb2*B1e.^2;

thetaEX(3) = thetaEX(2);
TR(3) = TR(2);
tm(3) = tm(2);

% Loop for each RF power and offset
for ii = 1:length([deltaCEST,deltaMT])
    for jj = 1:length(B1e)
                
        % Sat due to excitation
        C = diag([cos(thetaEX(jj)) 1]);
     
        % Calculate relaxation/exchange matrices
        kab = kba*M0(2)/M0(1);   
        if kba == 0
            kab = 0;
        end
        
        A = R1(1)*R1(2)+R1(1)*kba+R1(2)*kab;
        D = A+(R1(1)+kab)*Wb2(ii,jj)+(R1(2)+kba)*Wa(ii,jj)+...
            Wb2(ii,jj)*Wa(ii,jj);
        
        W = -diag([Wa(ii,jj) Wb2(ii,jj)]);
        
        % Steady State Magnetization
        Mss = D\[   M0(1)*(A+R1(1)*Wb2(ii,jj)); ...
                    M0(2)*(A+R1(2)*Wa(ii,jj))];
        
        % Relaxation Matrix 
        RL = [-R1(1)-kab            kba;...
                    kab         -R1(2)-kba];
        
        tr = TR(jj)-tm(jj)-ts;

        Es=expm(RL*ts); Er = expm(RL*tr);
        Em=expm((RL+W)*tm(jj));
        
        Mz(:,ii,jj) = (eye(2)-Es*Em*Er*C)\((Es*Em*(eye(2)-Er)+eye(2)-Es)...
            *M0+Es*(eye(2)-Em)*Mss);
        
    end
    
end


% Normalize according to same sequence w/out MT prepulse
Er0 = expm(RL*TR(jj));
Mz0 = (eye(3)-Er0*C)\(eye(3)-Er0)*M0;
Mzn = Mz ./ repmat(Mz0,[1 length([deltaCEST,deltaMT]) length(B1e)]);

% Extract free pool
Mz = squeeze(Mz(1,:,:))';
Mzn = squeeze(Mzn(1,:,:))';




