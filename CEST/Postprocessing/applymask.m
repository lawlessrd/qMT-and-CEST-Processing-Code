function [mm, ss, tt] = applymask(v, m)
% [mm,ss] = applymask(v,m)
%  v = 3-D volume
%  m = 3-D mask volume
%
%  mm = mean of masked volume
%  ss = STD of masked volume


mm = zeros(1, size(v,3));
ss = zeros(1, size(v,3));

for ii=1:size(v,3)
	vv = v(:,:,ii);
	mm(ii) = mean(vv(m==1));
	ss(ii) = std(vv(m==1));
    tt(ii) = sum(vv(m==1));
end
