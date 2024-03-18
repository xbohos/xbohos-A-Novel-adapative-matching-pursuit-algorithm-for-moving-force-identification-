function figure_plot_f
% Cite as: 'Xu B H, Yu L. A novel regularized adaptive matching pursuit for moving force identification using multiple criteria and prior knowledge[J]. 
% International Journal of Structural Stability and Dynamics, 2023, 23(10): 2350117.'
%20220517
close all
h=5e-3;%fs=200Hz
l=40;v=40;
n1=1;
h1=h/n1;
t=0:h1:l/v;ls=4;
[f11,f22]=ydhzdouble_force(t);
t0=0:h1:((l+ls)/v);
f1=zeros(length(t),1);f2=zeros(length(t),1);
f1(2:end-1)=f11(2:end-1);f2(2:end-1)=f22(2:end-1);
f1=[f1;zeros(length(t0)-length(f1),1)];f2=[zeros(length(t0)-length(f2),1);f2];
%%
load f_220512.mat
n_response=3;
f_OMP1=f_1(:,(n_response-1)*3+1);f_OMP1=[f_OMP1;zeros(length(t0)-length(f_OMP1),1)];
f_ROMP1=f_1(:,(n_response-1)*3+2);f_ROMP1=[f_ROMP1;zeros(length(t0)-length(f_ROMP1),1)];
f_huisuROMP1=f_1(:,(n_response-1)*3+3);f_huisuROMP1=[f_huisuROMP1;zeros(length(t0)-length(f_huisuROMP1),1)];
f_OMP2=f_2(:,(n_response-1)*3+1);f_OMP2=[zeros(length(t0)-length(f_OMP2)+1,1);f_OMP2(2:end)];
f_ROMP2=f_2(:,(n_response-1)*3+2);f_ROMP2=[zeros(length(t0)-length(f_ROMP2)+1,1);f_ROMP2(2:end)];
f_huisuROMP2=f_2(:,(n_response-1)*3+3);f_huisuROMP2=[zeros(length(t0)-length(f_huisuROMP2)+1,1);f_huisuROMP2(2:end)];
figure
hold on
plot(t0,f1,'Linewidth',2);plot(t0,f_OMP1,'--','Linewidth',1,'Color',[0.4660 0.6740 0.1880]);
plot(t0,f_ROMP1,':','Linewidth',2);plot(t0,f_huisuROMP1,'-.','Linewidth',2);
ylim([0 10]);
legend('True','OMP','ROMP','NRAMP');
%xlabel('Time/s');ylabel('Force/kN');
xlabel('Time/s');ylabel('Force/kN');
figure
hold on
plot(t0,f2,'Linewidth',2);plot(t0,f_OMP2,'--','Linewidth',1,'Color',[0.4660 0.6740 0.1880]);
plot(t0,f_ROMP2,':','Linewidth',2);plot(t0,f_huisuROMP2,'-.','Linewidth',2);
legend('True','OMP','ROMP','NRAMP');
ylim([0 25]);
xlabel('Time/s');ylabel(['Force/kN']);