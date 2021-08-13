function w1cw = w1cw_gauss ( degree, Tau, TR )
% this function calculates the cw w1 of a gauss pulse
% degree: mt pulse angle in degrees
% Tau : pulse duration
% TR: repetition time
% w1cw : in rad/s  
w1cw = zeros(length(degree),1);
for ii=1:length(degree)
    theta = degree(ii) * pi / 180;
    gam = 1; % then b1 = W1 in this function

    tstep = 1e-6; % 1ns step

    t= 0:tstep:Tau;

    T = sqrt(-(Tau/2)^2/log(0.01));
    B1a = theta/gam/Tau/0.4116;
    B1as = B1a*exp(-(t-Tau/2).^2./T^2);


    width = Tau * 42.67/256; % data from Varian pulse shape
    xc = Tau/2;
    B1b = theta/gam/Tau/0.4167; % calculated from Varian pulse        
    B1bs = B1b * exp ( -0.5* ((t- xc)/width).^2); % Varian Gaussian pulse shape

%     figure;
%     hold on;
%     plot( t, B1as, 'k-');
%     plot( t, B1bs, 'r-');
%     hold off;


    power = sum( B1bs .* B1bs) * tstep; % integral of power

    w1cw(ii)= sqrt ( power / TR);
end

end