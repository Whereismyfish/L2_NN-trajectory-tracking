


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
IC = [0 0 0 0.1 0 0 0 zeros(1,16) 0 0 0 0 0 0 0];
IC1 = [0 0 0 0.1 0 0 0 zeros(1,16) 0 0 0 0 0];
%=================

n1 = .5; n2 = .5;n3 = .5;

Error = odeset('RelTol',1e-3,'AbsTol',1e-4);
[T1,Y1] = ode45(@(t,y) NNoutput(t,y,m11,m22,m33,Xu,Yv,Nr,Xuu,Yvv,Nrr,k01,k02,k1,k2,k3,k4,k5,k6,delta,const_t,ganma,n1,n2,n3),Tspan,IC,Error);
[T2,Y2] = ode45(@(t,y) AUV_L2_noNN(t,y,m11,m22,m33,Xu,Yv,Nr,Xuu,Yvv,Nrr,k01,k02,k1,k2,k3,k4,k5,k6,delta,const_t,ganma,n1,n2,n3),Tspan,IC1,Error);


[~,tao_u,tao_r,e_u,e_r,dist1,dist3,x_d,y_d,Feu,Fer,x_e,y_e,psi_e,norm_Wu,norm_Wr,non1,non2] = cellfun(@(t,y) NNoutput(t,y.',m11,m22,m33,...
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

non1 = cell2mat(non1);non2 = cell2mat(non2);
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

figure(3)
subplot(211),plot(T1,non1-Feu,'linewidth',1.5),grid
ylabel('Surge motion','fontname','times','fontsize',14)
title('NN approximation errors','fontname','times','fontsize',14)

subplot(212),plot(T1,non2-Fer,'linewidth',1.5),grid
ylabel('Yaw motion','fontname','times','fontsize',14)
xlabel('\itt \rm(s)','fontname','times','fontsize',14),
figure(4)%%% modified on 2023-07-17
plot(T1,norm_Wu,T1,norm_Wr,'linewidth',1.5),grid
xlabel('\itt \rm(s)','fontname','times','fontsize',14),
ylabel('Norm of NN weights','fontname','times','fontsize',14)
legend('$\Vert \hat{W_u}\Vert$','$\Vert \hat{W_r}\Vert$','Interpreter','latex','fontname','Palatino Linotype','fontsize',12)



