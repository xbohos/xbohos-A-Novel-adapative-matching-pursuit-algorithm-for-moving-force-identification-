function main0515
% Cite as: Xu B H, Yu L. A novel regularized adaptive matching pursuit for moving force identification using multiple criteria and prior knowledge[J]. 
% International Journal of Structural Stability and Dynamics, 2023, 23(10): 2350117.
%2022.5.15
%% Measurement for SHM
close all
h=5e-3;%fs=200Hz
l=40;v=40;
n1=1;
h1=h/n1;
t=0:h1:l/v;
%load C_chongji.mat; 5%Noise
load C_chongji.mat;
load M_chongji_5Noise_PCD.mat;
load a_chongji_5Noise_PCD.mat;
[~,~,f11]=ydhzdouble_force(t);
M1=M21;M2=M22;M3=M23;a1=a21;a2=a22;a3=a23;
C4=C4./norm(M1);C5=C5./norm(M2);C6=C6./norm(M3);
C44=C44./norm(a1);C55=C55./norm(a2);C66=C66./norm(a3);
M1=M1./norm(M1);M2=M2./norm(M2);M3=M3./norm(M3);
a1=a1./norm(a1);a2=a2./norm(a2);a3=a3./norm(a3);
%Pre
fs=50;a=1e-10;alpha=0.5;
%the first article
%(1) 1/4m 1/2a
Phi=[C4;C55];
y=[M1;a2];
[error,f1,t11,a0,~]=ydhzsb_huisuROMP_first1(Phi,y,f11,fs,a,alpha);
error0(1:3,1)=error;
f_1(:,1:3)=f1;
t111(1:3,1)=t11';
static(1:3,1)=a0;
%%
%10%Noise
load C_chongji.mat;
load M_chongji_10Noise_PCD.mat;
load a_chongji_10Noise_PCD.mat;
M1=M21;M2=M22;M3=M23;a1=a21;a2=a22;a3=a23;
C4=C4./norm(M1);C5=C5./norm(M2);C6=C6./norm(M3);
C44=C44./norm(a1);C55=C55./norm(a2);C66=C66./norm(a3);
M1=M1./norm(M1);M2=M2./norm(M2);M3=M3./norm(M3);
a1=a1./norm(a1);a2=a2./norm(a2);a3=a3./norm(a3);
%Pre
fs=50;a=1e-10;alpha=0.5;
%the first article
%(1) 1/4m 1/2a
Phi=[C4;C55];
y=[M1;a2];
[error,f1,t11,a0,atom10]=ydhzsb_huisuROMP_first1(Phi,y,f11,fs,a,alpha);
error0(1:3,2)=error;
f_1(:,4:6)=f1;
t111(1:3,2)=t11';
static(1:3,2)=a0;
%%
%15%Noise
load C_chongji.mat;
load M_chongji_15Noise_PCD.mat;
load a_chongji_15Noise_PCD.mat;
M1=M21;M2=M22;M3=M23;a1=a21;a2=a22;a3=a23;
C4=C4./norm(M1);C5=C5./norm(M2);C6=C6./norm(M3);
C44=C44./norm(a1);C55=C55./norm(a2);C66=C66./norm(a3);
M1=M1./norm(M1);M2=M2./norm(M2);M3=M3./norm(M3);
a1=a1./norm(a1);a2=a2./norm(a2);a3=a3./norm(a3);
%Pre
fs=50;a=1e-10;alpha=0.5;
%the first article
%(1) 1/4m 1/2a
Phi=[C4;C55];
y=[M1;a2];
[error,f1,t11,a0,atom15]=ydhzsb_huisuROMP_first1(Phi,y,f11,fs,a,alpha);
error0(1:3,3)=error;
f_1(:,7:9)=f1;
t111(1:3,3)=t11';
static(1:3,3)=a0;
f1=zeros(length(f11),1);f1(2:end-1)=f11(2:end-1);
figure
hold on
plot(t,f1,'Linewidth',2);
plot(t,f_1(:,3),'--','Linewidth',2);
plot(t,f_1(:,6),':','Linewidth',2);
plot(t,f_1(:,9),'-.','Linewidth',2);
legend('True','Identified with 5% Noise','Identified with 10% Noise','Identified with 15% Noise');
xlabel('Time/s');ylabel('Force/kN');
end

function [error,f11,t11,a1,atom]=ydhzsb_huisuROMP_first1(Phi,y,f11,fs,a,alpha)
%21.5.21  %21.11.14 meeting  %22.5.3
%huisuROMP1 huisuROMP2
h=5e-3;%fs=200Hz
l=40;v=40;
n1=1;
h1=h/n1;
t=0:h1:l/v;
t1=t;nf1=length(t);nf2=length(t);
f1=zeros(length(t1),1);
f1(2:end-1)=f11(2:end-1);
%CS_OMP_origin
f1_OMP=zeros(length(t1),1);
[~,C2,~,~]=ydhzsb_redundantmatrix1(h,fs,nf1,nf2);
A1=Phi*C2;
K_OMP=round(size(A1,2)/2);a1=zeros(3,1);
tic
[x_OMP]=CS_OMP(y,A1,K_OMP);x_OMP1=C2*x_OMP;
f1_OMP(2:end-1)=x_OMP1;
t11(1)=toc;
a1(1)=x_OMP(1);
%
f1_ROMP=zeros(length(t1),1);
tic
x_ROMP=CS_ROMP(y,A1,K_OMP);x_ROMP1=C2*x_ROMP;
f1_ROMP(2:end-1)=x_ROMP1;
t11(2)=toc;
a1(2)=x_ROMP(1);
%CS_huisuRAMP
f1_huisuROMP2=zeros(nf1,1);
tic
x_I1=huisuRAMP(y,A1, a ,alpha);
f=C2*x_I1;
f1_huisuROMP2(2:end-1)=f;
t11(3)=toc;
f11=[f1_OMP f1_ROMP f1_huisuROMP2];
a1(3)=x_I1(1);
error=zeros(size(f11,2),1);
atom=[x_OMP x_ROMP x_I1];
for i=1:size(f11,2)
error(i,1)=ydhzsbwucha(f11(:,i),f1);
end
end

function error=ydhzsbwucha(fm,ft)
%2021.1.26 fm:fore measurement  ft:true force
error=norm(fm-ft)/norm(ft).*100;
format short
end

function Value_theta=huisuRAMP(y,A, a ,alpha)
%%220501 
%a : minmum of energy ;
%alpha : coefficient of force . 
[y_rows,y_columns] = size(y);
if y_rows<y_columns
    y = y';%y should be a column vector
end
[M,N] = size(A);
I = A(:,1);
%I = [];
position=[1];

A(:,1)=zeros(size(A,1),1);
Value_theta=zeros(N,round(size(A,2)/3));
Value_theta1=zeros(N,1);
a1= lsqr(I,y);
r_n=y-I*a1;
inter=1;
for ii=1:round(size(A,2)/3)
    u = A'*r_n;
    [~,pos] = selection1(u,a);
    if pos==0
        break
    end
    A_F = [I A(:,pos)];
    if size(A_F,2)>M
        
        break;
    end
    A(:,pos) = zeros(M,length(pos));
    pos=[position pos];pos_1=find(pos==1);x_0=zeros(length(pos),1);
    tol=1e-6;maxit=100;x_0(pos_1)=a1;
  theta_ls = lsqr(A_F,y,tol,maxit,[],[],x_0);
    [I1,Ii] = sort(abs(theta_ls),'descend');
    
    a1=theta_ls(pos_1);
    kk=1;
    pos1=[];
    for i1=1:length(Ii)
        if I1(i1)>alpha*a1
            if pos(Ii(i1))>1
                pos1(kk)=i1;
                kk=kk+1;
            end
        end
        if I1(i1)<alpha*a1
            break
        end
    end
    if kk>1
        Ii(pos1')=[];
    end
    I=A_F(:,Ii);
    x_I=theta_ls(Ii);
    position=[];
    position(1:size(x_I))=pos(Ii)';
    Value_theta(pos(Ii),inter)=x_I;
    Value_theta1(pos(Ii))=x_I;
    r_n = y - I*x_I;
    r(inter)=norm(r_n);
    inter=inter+1;
    if norm(r_n)<1e-6%Repeat the steps until r=0
        break;
    end
end
[~,pos]=sort(r);
Value_theta=Value_theta(:,pos(1));
end

function [val,pos] = selection1(u,a)
    productabs = abs(u);
    [productdes,indexproductdes] = sort(productabs,'descend');
    for ii = length(productdes):-1:1
        if productdes(ii)>1e-6
            break;
        end
    end
    Jval=productdes(1:ii);J=indexproductdes(1:ii);K=ii;
    %Identify:Choose a set J of the K biggest coordinates
    %Regularize:Among all subsets J0∈J with comparable coordinates
    MaxE = -1;
    for kk = 1:ii
        J0_tmp = zeros(1,K);iJ0 = 1;
        J0_tmp(iJ0) = J(kk);
        Energy = Jval(kk)^2;
        for mm = kk+1:K
            if Jval(kk)<2*Jval(mm)
                iJ0 = iJ0 + 1;
                J0_tmp(iJ0) = J(mm);
                Energy = Energy + Jval(mm)^2;
            else
                break;
            end
        end
        if Energy>MaxE
            J0 = J0_tmp(1:iJ0);
            MaxE = Energy;
        end
    end
    if MaxE<a
        disp('Error');
        pos=0;
        val=0;
    else
    pos =J0;
    val = productabs(J0);
    end
end