


close all
clc
clear

%%% Parameters %%%%%%
m11 = 215;%kg
m22 = 265;%kg
m33 = 80;%kg
Xu = 70;%kg/s
Yv = 100;%kg/s
Nr = 50;%kgm2/s
Xuu = 100;%kg/m
Yvv = 200;%kg/m
Nrr =100;%kgm2

k01 = .1; k02 = .1;
k1 = 20; k2 = 20;
k3 = 20; k4 = .001;
k5 = 20; k6 = .001;
delta = 5;
const_t = 0.1;
ganma = .1;

load matlab.mat

Tspan = out.t;
IC = [0 0 0 0.1 0 0 0 zeros(1,16)];
%=================

n1 = .5; n2 = .5;n3 = .5;
Error = odeset('RelTol',1e-3,'AbsTol',1e-4);
[T1,Y1] = ode45(@(t,y) AUV_NN_L2(t,y,m11,m22,m33,Xu,Yv,Nr,Xuu,Yvv,Nrr,k01,k02,k1,k2,k3,k4,k5,k6,delta,const_t,ganma,n1,n2,n3),Tspan,IC,Error);

[~,tao_u,tao_r,e_u,e_r,dist1,dist3,x_d,y_d,Feu,Fer,x_e,y_e,psi_e] = cellfun(@(t,y) AUV_NN_L2(t,y.',m11,m22,m33,...
    Xu,Yv,Nr,Xuu,Yvv,Nrr,k01,k02,k1,k2,k3,k4,k5,k6,delta,const_t,ganma,n1,n2,n3),num2cell(T1),num2cell(Y1,2),'uni',0);



%===== modified on 03-02-2017 =====
% figure(1)
% subplot(211),plot(T1,Y1(:,1),T1,cos(T1)),
% ylabel('\itx \rm(m)','fontname','times'),grid,title('x-coordinate'),
% 
% subplot(212),plot(T1,Y1(:,2),T1,sin(T1)),
% xlabel('\itt \rm(s)','fontname','times'),ylabel('\ity \rm(m)','fontname','times'),grid,title('y-coordinate')
x_d = cell2mat(x_d);y_d = cell2mat(y_d);
tao_u = cell2mat(tao_u);tao_r = cell2mat(tao_r);
Feu = cell2mat(Feu);Fer = cell2mat(Fer);
x_e = cell2mat(x_e);y_e = cell2mat(y_e);psi_e = cell2mat(psi_e);
e_u = cell2mat(e_u);e_r = cell2mat(e_r);

%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%
figure(1)
plot(Y1(:,1),Y1(:,2),out.x,out.y,x_d,y_d,'-.k','linewidth',1.5)
xlabel('\itx \rm(m)','fontname','times','fontsize',14)
ylabel('\ity \rm(m)','fontname','times','fontsize',14),grid
title('Trajectory tracking'),legend('Proposed','SMC','Desired')

figure(2)
subplot(311),plot(T1,x_e,out.t,out.xe,'linewidth',1.5),ylabel('\itx_e \rm(m)','fontname','times','fontsize',14),grid
legend('Proposed','SMC')
subplot(312),plot(T1,y_e,out.t,out.ye,'linewidth',1.5),ylabel('\ity_e \rm(m)','fontname','times','fontsize',14),grid
subplot(313),plot(T1,psi_e,out.t,out.psie,'linewidth',1.5),grid
xlabel('\itt \rm(s)','fontname','times','fontsize',14),
ylabel('\it\psi_e \rm(rad)','fontname','times','fontsize',14)

%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%
figure(3)
subplot(311),plot(T1,Y1(:,4),out.t,out.u,'linewidth',1.5),ylabel('\itu \rm(m/s)','fontname','times','fontsize',14),grid
legend('Proposed','SMC')
subplot(312),plot(T1,Y1(:,5),out.t,out.v,'linewidth',1.5),ylabel('\itv \rm(m/s)','fontname','times','fontsize',14),grid
subplot(313),plot(T1,Y1(:,6),out.t,out.r,'linewidth',1.5),grid
xlabel('\itt \rm(s)','fontname','times','fontsize',14),
ylabel('\itr \rm(rad/s)','fontname','times','fontsize',14)

% 
figure(4)
subplot(211),plot(T1,e_u,out.t,out.eu,'linewidth',1.5),ylabel('\ite_u \rm(m/s)','fontname','times','fontsize',14),grid
legend('Proposed','SMC')
subplot(212),plot(T1,e_r,out.t,out.er,'linewidth',1.5),grid
xlabel('\itt \rm(s)','fontname','times','fontsize',14),
ylabel('\ite_r \rm(rad/s)','fontname','times','fontsize',14)


% 
figure(5)
subplot(211),plot(T1,tao_u,out.t,out.tolu,'linewidth',1.5),ylim([-220 220]),grid
ylabel('\it\tau_u \rm(N)','fontname','times','fontsize',14)
legend('Proposed','SMC')


subplot(212),plot(T1,tao_r,out.t,out.tolr,'linewidth',1.5),ylim([-170 170]),grid
xlabel('\itt \rm(s)','fontname','times','fontsize',14),
ylabel('\it\tau_r \rm(Nm)','fontname','times','fontsize',14)

% figure(6)
% subplot(211),plot(T1,Feu),title('NN approximation in surge motion')
% subplot(212),plot(T1,Fer),title('NN approximation in yaw motion')
% % 
% 
% 
% 

