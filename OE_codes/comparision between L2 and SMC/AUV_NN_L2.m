
function [dy,tao_u,tao_r,e_u,e_r,dist1,dist3,x_d,y_d,Feu,Fer,x_e,y_e,psi_e] = AUV_NN_L2(t,y,m11,m22,m33,...
    Xu,Yv,Nr,Xuu,Yvv,Nrr,k01,k02,k1,k2,k3,k4,k5,k6,delta,const_t,ganma,n1,n2,n3)% modified on 03-02-2017
dy = zeros(23,1);

% y(1) ---- x coordiante in the earth-fixed frame
% y(2) ---- y coordiante in the earth-fixed frame
% y(3) ---- heading in the earth-fixed frame
% y(4) ---- surge speed in the body-fixed frame
% y(5) ---- sway speed in the body-fixed frame
% y(6) ---- yaw rate in the body-fixed frame

% y(7)----- lowpass filter output

% y(8) ---- first of disturbance observer state variables
% y(9) ---- second of disturbance observer state variables
% y(10) ---- third of disturbance observer state variables

%y(11)~y(18)--- NN updated weights

M_tans1 = [cos(y(3)) -sin(y(3)) 0
                sin(y(3))  cos(y(3)) 0
                0            0            1];
M_temp1 = M_tans1*[y(4);y(5);y(6)];

dy(1) = M_temp1(1);
dy(2) = M_temp1(2);
dy(3) = M_temp1(3);

u = y(4);v = y(5); r = y(6);

beta_drift = atan(v/(u));
 
% x_d = cos(t);y_d = sin(t);
% diff_x_d = -sin(t);diff_y_d = cos(t);

x_d = 0.5*t;y_d = 6 + 0.2*sin(0.2*t);
diff_x_d = 0.5;diff_y_d = 0.04*cos(0.2*t);
% 

d11 = Xu + Xuu*abs(u);
d22 = Yv + Yvv*abs(v);
d33 = Nr + Nrr*abs(r);
f_u = (m22*v*r - d11*u)/m11;
f_v = (-m11*u*r - d22*v)/m22;
f_r = ((m11-m22)*u*v - d33*r)/m33;


psi_F = atan(diff_y_d/diff_x_d);

Q = [cos(psi_F) sin(psi_F)
     -sin(psi_F) cos(psi_F)];
temp_xeye = Q*[(y(1)-x_d) (y(2)-y_d)]';

x_e = temp_xeye(1) ;
y_e = temp_xeye(2);

U_d = sqrt(diff_x_d^2 + diff_y_d^2);
psi_d = psi_F - beta_drift + atan(-y_e/delta);

psi_e = y(3) - psi_d;

dy(7) = (psi_d - y(7))/const_t;
d_psi_d = (psi_d - y(7))/const_t;

u_d = (-k01*x_e + U_d)*cos(beta_drift)*sqrt(delta^2 + y_e^2)/delta;
r_d = d_psi_d - k02*(y(3) - psi_d);

e_u = u - u_d;
e_r = r - r_d;

nninput = [u v r x_e y_e e_u e_r 1]'; % input of neural network

for i=1:length(nninput)
    fi(i)=1/(1+exp(-nninput(i)));
end

% for i=1:length(nninput) % leaky Relu function, on 2023-05-15
%     if nninput(i)>=0
%         fi(i)= nninput(i);
%     else
%         fi(i)= 0.001*nninput(i);
%     end
% end

% to describe the neural network weights
Wu = zeros(length(nninput),1);
Wr = zeros(length(nninput),1);
error = [e_u e_r]';

for j=1:length(nninput)
    dy(7+j) = k3*fi(j)*e_u - k4*norm(error)*y(7+j);
    Wu(j)=y(7+j);
end
for q=1:length(nninput)
    dy(7+length(Wu)+q) = k5*fi(q)*e_r - k6*norm(error)*y(7+length(Wu)+q);
    Wr(q)=y(7+length(Wu)+q);
end
tao_u = m11*(-Wu'*fi'- (1+1/4/ganma/ganma)*e_u - f_u - k1*e_u);
tao_r = m33*(-Wr'*fi'- (1+1/4/ganma/ganma)*e_r - f_r - k2*e_r);

Feu = Wu'*fi';
Fer = Wr'*fi';

% dist1 = n1*randn;
% dist2 = n2*randn;
% dist3 = n3*randn;

% dist1 = 12 + sin(0.8*t + pi/8) - 0.6*sin(0.5*pi*t);
% dist2 = n2*randn;
% dist3 = 5 + sin(0.8*t + pi/6) - 0.5*sin(0.3*pi*t);

dist1 = 0.3*sin(0.2*t);
dist2 = 0.1*sin(0.2*t);
dist3 = 0.1*cos(0.3*t);



B_tao_u  = 200;
if tao_u>B_tao_u
    tao_u = B_tao_u;
elseif tao_u<-B_tao_u
    tao_u = -B_tao_u;
elseif tao_u>=-B_tao_u && tao_u<=B_tao_u
    tao_u = B_tao_u*tanh(tao_u/B_tao_u);
end

% B_tao_u  = 200;
% if tao_u>B_tao_u
%     tao_u = B_tao_u;
% elseif tao_u<-B_tao_u
%     tao_u = -B_tao_u;
% end


B_tao_r  = 150;
if tao_r>B_tao_r
    tao_r = B_tao_r;
elseif tao_r<-B_tao_r
    tao_r = -B_tao_r;
elseif tao_r>=-B_tao_r && tao_r<=B_tao_r
    tao_r = B_tao_r*tanh(tao_r/B_tao_r);
end

% B_tao_r  = 150;
% if tao_r>B_tao_r
%     tao_r = B_tao_r;
% elseif tao_r<-B_tao_r
%     tao_r = -B_tao_r;
% end


%============  modified on 03-02-2017 ======

%============  modified on 03-02-2017 ======

dy(4) = f_u + tao_u/m11 + dist1;
dy(5) = f_v + dist2;
dy(6) = f_r + tao_r/m33 + dist3;



        

















