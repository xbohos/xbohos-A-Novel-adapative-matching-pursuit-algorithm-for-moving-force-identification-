function [C3,C5,n1,n2]=ydhzsb_redundantmatrix1(h,fs,nf1,nf2)
%21.3.21
% Ref: [1] C. D. Pan et al., Moving force identi¯cation based on redundant concatenated%
%dictionary and weighted l1-norm regularization, Mech. Syst. Signal Process. 98 (2018)
%32–49.
% [2] Xu B H, Yu L. A novel regularized adaptive matching pursuit for moving force identification using multiple criteria and prior knowledge[J]. 
% International Journal of Structural Stability and Dynamics, 2023, 23(10): 2350117.
m1=nf1;m2=nf2;t1=0:1/(m1-1):1;t2=0:1/(m2-1):1;
T1=h*nf1;T2=h*nf2;n1=round(2*T1*fs+0.5);n2=round(2*T2*fs+0.5);
C_11_s=zeros(m1,n1);C_22_s=zeros(m2,n2);C_11_r=zeros(m1,n1);C_22_r=zeros(m2,n2);C_11_c=zeros(m1,n1+1);C_22_c=zeros(m2,n2+1);
C_11_c(:,1)=(2*T1)^0.5/T1.*ones(m1,1);
for i=2:n1+1
    for j=1:m1
        C_11_c(j,i)=(2*T1)^0.5/T1*cos((i-1)*pi*t1(j));
     
    end
end
for i=1:n1
    for j=1:m1
      C_11_s(j,i)=(2*T1)^0.5/T1*sin(i*pi*t1(j));
    
    end
end
C_22_c(:,1)=(T2)^0.5/T2.*ones(m2,1);
%C_22_c(:,1)=(nf2*T2)^0.5/T2.*ones(m2,1);
for i=2:n2+1
    for j=1:m2
      C_22_c(j,i)=(2*T2)^0.5/T2*cos((i-1)*pi*t2(j));
    end
end
for i=1:n2
    for j=1:m2
     C_22_s(j,i)=(2*T2)^0.5/T2*sin(i*pi*t2(j));
    end
end

for i=1:n1
    for j=1:m1
        if ((j-1)<i||(j-1)==i) && i<j+1
        C_11_r(j,i)=((nf1*T1)^0.5)/T1;
        end
    end
end
for i=1:n2
    for j=1:m2
        if ((j-1)<i||(j-1)==i)==i&&i<j+1
        C_22_r(j,i)=((nf2*T2)^0.5)/T2;
        end
    end
end
C_11_s=C_11_s(2:end-1,:);C_11_r=C_11_r(2:end-1,:);C_11_c=C_11_c(2:end-1,:);
C_22_s=C_22_s(2:end-1,:);C_22_r=C_22_r(2:end-1,:);C_22_c=C_22_c(2:end-1,:);
C3=zeros(m1+m2-4,n1*3+n2*3+2);
C3(1:m1-2,1:3*n1+1)=[ C_11_c C_11_s C_11_r];
C3(m1-1:m1+m2-4,3*n1+2:end)=[C_22_c C_22_s C_22_r];
C5=[ C_11_c C_11_s C_11_r];
