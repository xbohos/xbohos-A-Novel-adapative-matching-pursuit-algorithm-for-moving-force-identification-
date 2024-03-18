function main0514_1
%2022.5.14
% Cite as: Xu B H, Yu L. A novel regularized adaptive matching pursuit for moving force identification using multiple criteria and prior knowledge[J]. 
% International Journal of Structural Stability and Dynamics, 2023, 23(10): 2350117.'
% Validate the huisuRAMP in numerical simualtion
%3.7   %4.2  
%
close all
h=5e-3;%fs=200Hz
l=40;v=40;
n1=1;
h1=h/n1;
t=0:h1:l/v;
load C.mat;
load a_sin0520_10Noise_PCD.mat;
load M_sin0520_10Noise_PCD.mat;
[f11,f22,~]=ydhzdouble_force(t);
M1=M21;M2=M22;M3=M23;a1=a21;a2=a22;a3=a23;
C4=C4./norm(M1);C5=C5./norm(M2);C6=C6./norm(M3);
C44=C44./norm(a1);C55=C55./norm(a2);C66=C66./norm(a3);
M1=M1./norm(M1);M2=M2./norm(M2);M3=M3./norm(M3);
a1=a1./norm(a1);a2=a2./norm(a2);a3=a3./norm(a3);
%Pre
fs=50;a=1e-10;alpha=0.2;
subsystems = {
    {C4, C55},                          
    {C5, C55},                          
    {C4, C55, C66},                    
    {C4, C5, C6},                       
    {C4, C66, C5},                     
    {C44, C5, C66},                     
    {C4, C6, C44, C66},               
    {C4, C44, C55, C66}                 
};

y_subsystems = {
    {M1, a2},                         
    {M2, a2},                           
    {M1, a2, a3},                      
    {M1, M2, M3},                       
    {M1, a3, M2},                       
    {a1, M2, a3},                      
    {M1, M3, a1, a3},                   
    {M1, a1, a2, a3}                   
};


for i = 1:numel(subsystems)
    Phi = vertcat(subsystems{i}{:});
    y = vertcat(y_subsystems{i}{:});
    [error, f1, f2, t11] = ydhzsb_huisuROMP_first1(Phi, y, f11, f22, fs, a, alpha);
    error0(1:3,(i - 1) * 2 + 1:i * 2) = error;
    f_1(:,(i - 1) * 3 + 1:i * 3) = f1;
    f_2(:,(i - 1) * 3 + 1:i * 3) = f2;
    t111((i - 1) * 3 + 1:i * 3) = t11';
end

save('error_220512.mat','error0');
save('f_220512.mat','f_1','f_2');
save('t_220512.mat','t111');
end

function [error,f11,f22,t11]=ydhzsb_huisuROMP_first1(Phi,y,f11,f22,fs,a,alpha)
%21.5.21  %21.11.14 meeting  %22.5.3
%huisuROMP1 huisuROMP2
h=5e-3;%fs=200Hz
l=40;v=40;
n1=1;
h1=h/n1;
t=0:h1:l/v;
t1=t;t2=t;
f1=zeros(length(t1),1);
f1(2:end-1)=f11(2:end-1);
f2=zeros(length(t2),1);
f2(2:end-1)=f22(2:end-1);nf1=length(t);nf2=length(t);
%CS_OMP_origin
tic
f1_OMP_origin=zeros(length(t1),1);f2_OMP_origin=zeros(length(t2),1);
[C2,~,~,~]=ydhzsb_redundantmatrix1(h,fs,nf1,nf2);
A1=Phi*C2;
K_OMP=size(A1,2)/2;
[x_OMP]=CS_OMP(y,A1,K_OMP);x_OMP1=C2*x_OMP;
f1_OMP_origin(2:end-1)=x_OMP1(1:size(Phi,2)/2);f2_OMP_origin(2:end-1)=x_OMP1(size(Phi,2)/2+1:end);
t11(1)=toc;
%CS_ROMP
tic
f1_ROMP_origin=zeros(length(t1),1);f2_ROMP_origin=zeros(length(t2),1);
[C2,~,~,~]=ydhzsb_redundantmatrix1(h,fs,nf1,nf2);
K_OMP=size(A1,2)/2;
[x_ROMP]=CS_ROMP(y,A1,K_OMP);x_ROMP1=C2*x_ROMP;
f1_ROMP_origin(2:end-1)=x_ROMP1(1:size(Phi,2)/2);f2_ROMP_origin(2:end-1)=x_ROMP1(size(Phi,2)/2+1:end);
t11(2)=toc;
%CS_huisuRAMP
f1_huisuROMP2=zeros(nf1,1);f2_huisuROMP2=zeros(nf2,1);
tic
x_I1=huisuRAMP(y,A1, a ,alpha );
f=C2*x_I1;
f1_huisuROMP2(2:end-1)=f(1:nf1-2);f2_huisuROMP2(2:end-1)=f(nf1-1:end);
t11(3)=toc;
f11=[f1_OMP_origin f1_ROMP_origin f1_huisuROMP2];
f22=[f2_OMP_origin f2_ROMP_origin f2_huisuROMP2];
error=zeros(size(f11,2),2);
for i=1:size(f11,2)
error(i,1)=ydhzsbwucha(f11(:,i),f1);
error(i,2)=ydhzsbwucha(f22(:,i),f2);
end

end

function error=ydhzsbwucha(fm,ft)
%2021.1.26 fm:fore measurement  ft:true force
error=norm(fm-ft)/norm(ft).*100;
format short
end




function Value_theta=huisuRAMP(y,A, a ,alpha )
%%220501 
%a : minmum of energy ;
%alpha : coefficient of force . 
[y_rows,y_columns] = size(y);
if y_rows<y_columns
    y = y';%y should be a column vector
end
[M,N] = size(A);
I = [A(:,1) A(:,size(A,2)/2+1)];
position=[1 size(A,2)/2+1];
A(:,1)=zeros(size(A,1),1);A(:,size(A,2)/2+1)=zeros(size(A,1),1);
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
        break
    end
    A(:,pos) = zeros(M,length(pos));
    %%Consider Biconjugate gradient stabilized method
   
    pos=[position pos];pos_1=find(pos==1);pos_2=find(pos==size(A,2)/2+1);x_0=zeros(length(pos),1);
    tol=1e-6;maxit=100;x_0(pos_1)=a1(1);x_0(pos_2)=a1(2);
  theta_ls = lsqr(A_F,y,tol,maxit,[],[],x_0);
    [I1,Ii] = sort(abs(theta_ls),'descend');
    
    a1(1)=theta_ls(pos_1);a1(2)=theta_ls(pos_2);
    kk=1;
    pos1=[];
    for i1=1:length(Ii)
        if I1(i1)>alpha*a1(2)
            if pos(Ii(i1))>(size(A,2)/2+1)
                pos1(kk)=i1;
                kk=kk+1;
            end
        end
        if I1(i1)>alpha*a1(1)
           if pos(Ii(i1))<(size(A,2)/2+1)&&pos(Ii(i1))>1
                pos1(kk)=i1;
                kk=kk+1;
            end
        end
        if I1(i1)<alpha*a1(1)
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
        break
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
    %Regularize:Among all subsets J0¡ÊJ with comparable coordinates
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