function [f1,f2,fchongji]=ydhzdouble_force(t)
% Ref: [1] C. D. Pan et al., Moving force identi¯cation based on redundant concatenated%
%dictionary and weighted l1-norm regularization, Mech. Syst. Signal Process. 98 (2018)
%32–49.
% [2] Xu B H, Yu L. A novel regularized adaptive matching pursuit for moving force identification using multiple criteria and prior knowledge[J]. 
% International Journal of Structural Stability and Dynamics, 2023, 23(10): 2350117.
f1=5.*(1+0.1*sin(10*pi*t)+0.05*sin(40*pi*t));
f2=20.*(1-0.1*sin(10*pi*t)+0.05*sin(50*pi*t));
fchongji=zeros(length(t),1);
for i=1:length(t)
    if t(i)<0.6
     fchongji(i)=40*(1+0.3*sin(25*pi*t(i))+0.2*sin(60*pi*t(i)));
    else      
fchongji(i)=40*(1+0.3*sin(25*pi*t(i))+0.2*sin(60*pi*t(i))+3*exp(-35*(t(i)-0.6)).*sin(125*(t(i)-0.6)));
    end
end