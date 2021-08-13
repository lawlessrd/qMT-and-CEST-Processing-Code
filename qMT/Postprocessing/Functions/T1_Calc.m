function [T1obs] = T1_Calc(T1flip,MFA,MFA_Data,corrB1,T1TR)

del_flip = T1flip/MFA;
thetaT1 = T1flip:-del_flip:del_flip; % deg

yval = MFA_Data./sind(thetaT1*corrB1);
xval = MFA_Data./tand(thetaT1*corrB1);
pfit = polyfit(xval,yval,1);

E1 = pfit(1);

T1obs = -T1TR/log(E1);
end

