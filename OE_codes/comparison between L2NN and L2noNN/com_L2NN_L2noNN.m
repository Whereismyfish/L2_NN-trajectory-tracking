


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

L0 = diag([10 10 10]);

Tspan = 0:0.1:200;
IC = [0 0 0 0.1 0 0 0 zeros(1,16) 0 0 0 0 0];
IC1 = [0 0 0 0.1 0 0 0 zeros(1,16) 0 0 0 0 0 0 0 0];
%=================

n1 = .5; n2 = .5;n3 = .5;

Error = odeset('RelTol',1e-3,'AbsTol',1e-4);
[T1,Y1] = ode45(@(t,y) AUV_L2NN(t,y,m11,m22,m33,Xu,Yv,Nr,Xuu,Yvv,Nrr,k01,k02,k1,k2,k3,k4,k5,k6,delta,const_t,ganma,n1,n2,n3),Tspan,IC,Error);
[T2,Y2] = ode45(@(t,y) AUV_L2_noNN(t,y,m11,m22,m33,Xu,Yv,Nr,Xuu,Yvv,Nrr,k01,k02,k1,k2,k3,k4,k5,k6,delta,const_t,ganma,n1,n2,n3),Tspan,IC,Error);
kb = 7;
% [T4,Y4] = ode45(@(t,y) AUV_NN_L2_BLF(t,y,m11,m22,m33,Xu,Yv,Nr,Xuu,Yvv,Nrr,k01,k02,k1,k2,k3,k4,k5,k6,delta,const_t,ganma,n1,n2,n3,kb),Tspan,IC,Error);


[~,tao_u,tao_r,e_u,e_r,dist1,dist3,x_d,y_d,Feu,Fer,x_e,y_e,psi_e,norm_Wu,norm_Wr] = cellfun(@(t,y) AUV_L2NN(t,y.',m11,m22,m33,...
    Xu,Yv,Nr,Xuu,Yvv,Nrr,k01,k02,k1,k2,k3,k4,k5,k6,delta,const_t,ganma,n1,n2,n3),num2cell(T1),num2cell(Y1,2),'uni',0);

[~,tao_u1,tao_r1,e_u1,e_r1,dist11,dist31,x_d1,y_d1,Feu1,Fer1,x_e1,y_e1,psi_e1] = cellfun(@(t,y) AUV_L2_noNN(t,y.',m11,m22,m33,...
    Xu,Yv,Nr,Xuu,Yvv,Nrr,k01,k02,k1,k2,k3,k4,k5,k6,delta,const_t,ganma,n1,n2,n3),num2cell(T2),num2cell(Y2,2),'uni',0);


%===== modified on 03-02-2017 =====
% figure(1)
% subplot(211),plot(T1,Y1(:,1),T1,cos(T1)),
% ylabel('\itx \rm(m)','fontname','times'),grid,title('x-coordinate'),
% 
% subplot(212),plot(T1,Y1(:,2),T1,sin(T1)),
% xlabel('\itt \rm(s)','fontname','times'),ylabel('\ity \rm(m)','fontname','times'),grid,title('y-coordinate')
x_d = cell2mat(x_d);y_d = cell2mat(y_d);
x_d1 = cell2mat(x_d1);y_d1 = cell2mat(y_d1);

tao_u = cell2mat(tao_u);tao_r = cell2mat(tao_r);
tao_u1 = cell2mat(tao_u1);tao_r1 = cell2mat(tao_r1);

Feu = cell2mat(Feu);Fer = cell2mat(Fer);

x_e = cell2mat(x_e);y_e = cell2mat(y_e);psi_e = cell2mat(psi_e);
x_e1 = cell2mat(x_e1);y_e1 = cell2mat(y_e1);psi_e1 = cell2mat(psi_e1);

e_u = cell2mat(e_u);e_r = cell2mat(e_r);
e_u1 = cell2mat(e_u1);e_r1 = cell2mat(e_r1);

dist1 = cell2mat(dist1);dist3 = cell2mat(dist3);
dist11 = cell2mat(dist11);dist31 = cell2mat(dist31);

norm_Wu = cell2mat(norm_Wu);norm_Wr = cell2mat(norm_Wr);

figure(1)
plot(Y1(:,1),Y1(:,2),'m',Y2(:,1),Y2(:,2),'b',x_d,y_d,'k-.','linewidth',1.5),hold on
xlabel('\itx \rm(m)','fontname','times','fontsize',14)
ylabel('\ity \rm(m)','fontname','times','fontsize',14),grid
title('Trajectory tracking'),legend('L2 and NN based','L2 without NN based','Desired')

figure(2)
subplot(311),plot(T1,x_e,'m',T2,x_e1,'b','linewidth',1.5),ylabel('\itx_e \rm(m)','fontname','times','fontsize',14),grid
legend('L2 and NN based','L2 without NN based')
subplot(312),plot(T1,y_e,'m',T2,y_e1,'b','linewidth',1.5),ylabel('\ity_e \rm(m)','fontname','times','fontsize',14),grid
subplot(313),plot(T1,psi_e,'m',T2,psi_e1,'b','linewidth',1.5),grid
xlabel('\itt \rm(s)','fontname','times','fontsize',14),
ylabel('\it\psi_e \rm(rad)','fontname','times','fontsize',14)

figure(3)%%% modified on 2023-07-17
plot(T1,norm_Wu,T1,norm_Wr,'linewidth',1.5),grid
figure(7)
subplot(221),plot(T1(1:100),psi_e(1:100),'m',T2(1:100),psi_e1(1:100),'b','linewidth',1.5),grid
xlabel('\itt \rm(s)','fontname','times','fontsize',14),
ylabel('\it\psi_e \rm(rad)','fontname','times','fontsize',14)

subplot(222),plot(T1(1000:1200),psi_e(1000:1200),'m',T2(1000:1200),psi_e1(1000:1200),'b','linewidth',1.5),grid
xlabel('\itt \rm(s)','fontname','times','fontsize',14),
ylabel('\it\psi_e \rm(rad)','fontname','times','fontsize',14)

subplot(212),plot(T1(1300:2000),psi_e(1300:2000),'m',T2(1300:2000),psi_e1(1300:2000),'b','linewidth',1.5),grid
ylabel('\it\psi_e \rm(rad)','fontname','times','fontsize',14)
xlabel('\itt \rm(s)','fontname','times','fontsize',14),




figure(3)
subplot(311),plot(T1,Y1(:,4),'m',T2,Y2(:,4),'b','linewidth',1.5),ylabel('\itu \rm(m/s)','fontname','times','fontsize',14),grid
legend('Proposed','Do-based method')
subplot(312),plot(T1,Y1(:,5),'m',T2,Y2(:,5),'b','linewidth',1.5),ylabel('\itv \rm(m/s)','fontname','times','fontsize',14),grid
subplot(313),plot(T1,Y1(:,6),'m',T2,Y2(:,6),'b','linewidth',1.5),xlabel('\itt \rm(s)','fontname','times','fontsize',14),grid
ylabel('\itr \rm(rad/s)','fontname','times','fontsize',14)
% 
figure(8)
subplot(221),plot(T1(1:100),Y1(1:100,4),'m',T2(1:100),Y2(1:100,4),'b','linewidth',1.5),grid
xlabel('\itt \rm(s)','fontname','times','fontsize',14),
ylabel('\itu \rm(m/s)','fontname','times','fontsize',14)

subplot(222),plot(T1(1000:1200),Y1(1000:1200,4),'m',T2(1000:1200),Y2(1000:1200,4),'b','linewidth',1.5),grid
xlabel('\itt \rm(s)','fontname','times','fontsize',14),
ylabel('\itu \rm(m/s)','fontname','times','fontsize',14)

subplot(212),plot(T1(1300:2000),Y1(1300:2000,4),'m',T2(1300:2000),Y2(1300:2000,4),'b','linewidth',1.5),grid
ylabel('\itu \rm(m/s)','fontname','times','fontsize',14)
xlabel('\itt \rm(s)','fontname','times','fontsize',14),

figure(9)%%% small plot of v velocity
subplot(221),plot(T1(1:100),Y1(1:100,5),'m',T2(1:100),Y2(1:100,5),'b','linewidth',1.5),grid
xlabel('\itt \rm(s)','fontname','times','fontsize',14),
ylabel('\itv \rm(m/s)','fontname','times','fontsize',14)

subplot(222),plot(T1(1000:1200),Y1(1000:1200,5),'m',T2(1000:1200),Y2(1000:1200,5),'b','linewidth',1.5),grid
xlabel('\itt \rm(s)','fontname','times','fontsize',14),
ylabel('\itv \rm(m/s)','fontname','times','fontsize',14)

subplot(212),plot(T1(1300:2000),Y1(1300:2000,5),'m',T2(1300:2000),Y2(1300:2000,5),'b','linewidth',1.5),grid
ylabel('\itv \rm(m/s)','fontname','times','fontsize',14)
xlabel('\itt \rm(s)','fontname','times','fontsize',14),

figure(10)%%% small plot of r velocity
subplot(221),plot(T1(1:100),Y1(1:100,6),'m',T2(1:100),Y2(1:100,6),'b','linewidth',1.5),grid
xlabel('\itt \rm(s)','fontname','times','fontsize',14),
ylabel('\itr \rm(rad/s)','fontname','times','fontsize',14)

subplot(222),plot(T1(1000:1200),Y1(1000:1200,6),'m',T2(1000:1200),Y2(1000:1200,6),'b','linewidth',1.5),grid
xlabel('\itt \rm(s)','fontname','times','fontsize',14),
ylabel('\itr \rm(rad/s)','fontname','times','fontsize',14)

subplot(212),plot(T1(1300:2000),Y1(1300:2000,6),'m',T2(1300:2000),Y2(1300:2000,6),'b','linewidth',1.5),grid
ylabel('\itr \rm(rad/s)','fontname','times','fontsize',14)
xlabel('\itt \rm(s)','fontname','times','fontsize',14),

figure(4)
subplot(211),plot(T1,e_u,'m',T2,e_u1,'b','linewidth',1.5),ylabel('\ite_u \rm(m/s)','fontname','times','fontsize',14),grid
legend('Proposed','Do-based method')
subplot(212),plot(T1,e_r,'m',T2,e_r1,'b','linewidth',1.5),grid
xlabel('\itt \rm(s)','fontname','times','fontsize',14),
ylabel('\ite_r \rm(rad/s)','fontname','times','fontsize',14)

figure(13)%%% small plot of velocity tracking errors
subplot(221),plot(T1(1:100),e_u(1:100),'m',T2(1:100),e_u1(1:100),'b','linewidth',1.5),grid
xlabel('\itt \rm(s)','fontname','times','fontsize',14),
ylabel('\ite_u \rm(m/s)','fontname','times','fontsize',14)

subplot(222),plot(T1(1000:1200),e_u(1000:1200),'m',T2(1000:1200),e_u1(1000:1200),'b','linewidth',1.5),grid
xlabel('\itt \rm(s)','fontname','times','fontsize',14),
ylabel('\ite_u \rm(m/s)','fontname','times','fontsize',14)

subplot(212),plot(T1(1300:2000),e_u(1300:2000),'m',T2(1300:2000),e_u1(1300:2000),'b','linewidth',1.5),grid
ylabel('\ite_u \rm(m/s)','fontname','times','fontsize',14)
xlabel('\itt \rm(s)','fontname','times','fontsize',14),

figure(14)%%% small plot of velocity tracking errors
subplot(221),plot(T1(1:100),e_r(1:100),'m',T2(1:100),e_r1(1:100),'b','linewidth',1.5),grid
xlabel('\itt \rm(s)','fontname','times','fontsize',14),
ylabel('\ite_r \rm(rad/s)','fontname','times','fontsize',14)

subplot(222),plot(T1(1000:1200),e_r(1000:1200),'m',T2(1000:1200),e_r1(1000:1200),'b','linewidth',1.5),grid
xlabel('\itt \rm(s)','fontname','times','fontsize',14),
ylabel('\ite_r \rm(rad/s)','fontname','times','fontsize',14)

subplot(212),plot(T1(1300:2000),e_r(1300:2000),'m',T2(1300:2000),e_r1(1300:2000),'b','linewidth',1.5),grid
ylabel('\ite_r \rm(rad/s)','fontname','times','fontsize',14)
xlabel('\itt \rm(s)','fontname','times','fontsize',14),

figure(5)
subplot(211),plot(T1,tao_u,'m',T2,tao_u1,'b','linewidth',1.5),ylabel('\it\tau_u \rm(N)','fontname','times','fontsize',14),ylim([-50 210]),grid
legend('Proposed','Do-based method')
subplot(212),plot(T1,tao_r,'m',T2,tao_r1,'b','linewidth',1.5),ylim([-160 160]),grid
xlabel('\itt \rm(s)','fontname','times','fontsize',14),
ylabel('\it\tau_r \rm(Nm)','fontname','times','fontsize',14),hold on

figure(11)%%% small plot of control inputs
subplot(221),plot(T1(1:100),tao_u(1:100),'m',T2(1:100),tao_u1(1:100),'b','linewidth',1.5),grid
xlabel('\itt \rm(s)','fontname','times','fontsize',14),
ylabel('\it\tau_u \rm(N)','fontname','times','fontsize',14)

subplot(222),plot(T1(1000:1200),tao_u(1000:1200),'m',T2(1000:1200),tao_u1(1000:1200),'b','linewidth',1.5),grid
xlabel('\itt \rm(s)','fontname','times','fontsize',14),
ylabel('\it\tau_u \rm(N)','fontname','times','fontsize',14)

subplot(212),plot(T1(1300:2000),tao_u(1300:2000),'m',T2(1300:2000),tao_u1(1300:2000),'b','linewidth',1.5),grid
ylabel('\it\tau_u \rm(N)','fontname','times','fontsize',14)
xlabel('\itt \rm(s)','fontname','times','fontsize',14),

figure(12)%%% small plot of control inputs
subplot(221),plot(T1(1:100),tao_r(1:100),'m',T2(1:100),tao_r1(1:100),'b','linewidth',1.5),grid
xlabel('\itt \rm(s)','fontname','times','fontsize',14),
ylabel('\it\tau_r \rm(Nm)','fontname','times','fontsize',14)

subplot(222),plot(T1(1000:1200),tao_r(1000:1200),'m',T2(1000:1200),tao_r1(1000:1200),'b','linewidth',1.5),grid
xlabel('\itt \rm(s)','fontname','times','fontsize',14),
ylabel('\it\tau_r \rm(Nm)','fontname','times','fontsize',14)

subplot(212),plot(T1(1300:2000),tao_r(1300:2000),'m',T2(1300:2000),tao_r1(1300:2000),'b','linewidth',1.5),grid
ylabel('\it\tau_r \rm(Nm)','fontname','times','fontsize',14)
xlabel('\itt \rm(s)','fontname','times','fontsize',14),

% figure(6)
% subplot(311),plot(T2,dist11,'linewidth',1.5),title('Disturbances')
% ylabel('\itd_1','fontname','times','fontsize',14),grid
% % legend('\itd_1','Estimate of \itd_1')
% subplot(312),plot(T2,dist21,'linewidth',1.5),ylabel('\itd_2','fontname','times','fontsize',14),grid
% subplot(313),plot(T2,dist31,'linewidth',1.5),
% xlabel('\itt \rm(s)','fontname','times','fontsize',14),
% ylabel('\itd_3','fontname','times','fontsize',14),grid

% 
% figure(6)
% subplot(211),plot(T1,Feu),title('NN approximation in surge motion')
% subplot(212),plot(T1,Fer),title('NN approximation in yaw motion')
% 
%%% position error analysis %%%%%


