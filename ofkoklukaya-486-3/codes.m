%% initial

clc
clear
close all
load hw2_data.mat;

O1=zeros(3,25000);
O2=zeros(3,25000);
O3=zeros(3,25000);
sun_orb=zeros(3,25000);
sun_body=zeros(3,25000);
A_oi=zeros(3,3,25000);
A_bi=zeros(3,3,25000);
q_bi=zeros(4,25000);
wdot_bi=zeros(3,25000);
w_bi=zeros(3,25000);
w_i=zeros(3,25000);
origin=zeros(3,25000);
N_0=zeros(3,1);
t=zeros(1,25000);
rpm2rads=2*pi/60;
w_bi(:,1)=[0 0 60*rpm2rads];
I3=eye(3);


sun_body(:,1)=[0; 0; 1];


% A_bi(:,:,1)=sun_body(:,1)*transpose(sun_eci(:,1)); %wrong approach

A_bi(:,:,1)=[0.0464 0.999 0;...
    -0.02 0.001 -1;...
    -0.999 -0.046 0.02];
 q_bi(:,1)=dcm2quat(A_bi(:,:,1));
  q_bi(:,1)=[q_bi(2,1),q_bi(3,1),q_bi(4,1),q_bi(1,1)]; %changing to our notation


%% Q1

dt=0.01;
t(1)=1;

limit=1000;


J=[0.018,0,0;...
    0,0.018,0;...
    0,0,0.018];

for s=1:limit

    q_bi(:,s)=[q_bi(4,s),q_bi(1,s),q_bi(2,s),q_bi(3,s)];%
    A_bi(:,:,s)=quat2dcm(transpose(q_bi(:,s)));
    q_bi(:,s)=[q_bi(2,s),q_bi(3,s),q_bi(4,s),q_bi(1,s)];%


    qk1=0.5*omgof(w_bi(:,s))*q_bi(:,s); %rk4
    wk1=inv(J)*(N_0-(cross(w_bi(:,s),(J*w_bi(:,s)))));
    qrk=q_bi(:,s)+qk1*(dt/2);
    wrk=w_bi(:,s)+wk1*(dt/2);
    qk2=0.5*omgof(wrk)*qrk;
    wk2=inv(J)*(N_0-(cross(wrk,(J*wrk))));
    qrk=q_bi(:,s)+qk2*(dt/2);
    wrk=w_bi(:,s)+wk2*(dt/2);
    qk3=0.5*omgof(wrk)*qrk;
    wk3=inv(J)*(N_0-(cross(wrk,(J*wrk))));
    qrk=q_bi(:,s)+qk2*(dt);
    wrk=w_bi(:,s)+wk2*(dt);
    qk4=0.5*omgof(wrk)*qrk;
    wk4=inv(J)*(N_0-(cross(wrk,(J*wrk))));
    q_bi(:,s+1)=q_bi(:,s)+((qk1+2*qk2+2*qk3+qk4)/6)*dt;
    w_bi(:,s+1)=w_bi(:,s)+((wk1+2*wk2+2*wk3+wk4)/6)*dt;

    q_bi(:,s+1)=q_bi(:,s+1)/norm(q_bi(:,s+1));
    t(s)=(s-1)*dt;
end

for i=1:limit
    w_i(:,i)=transpose(A_bi(:,:,i))*w_bi(:,i);
    spinaxis(:,i)=w_i(:,i)/norm(w_i(:,i));
end


figure(1);
subplot(3,2,1)
plot(t(1:limit),spinaxis(1,1:limit));
title("x-components of the spin axis in the inertial frame");
xlabel("t(s)");
ylim([-1 1])
subplot(3,2,3)
plot(t(1:limit),spinaxis(2,1:limit));
title("y-components of the spin axis in the inertial frame");
xlabel("t(s)");
ylim([-1 1])
subplot(3,2,5)
plot(t(1:limit),spinaxis(3,1:limit));
title("z-components of the spin axis in the inertial frame");
xlabel("t(s)");
ylim([-1 1])
hold on;
subplot(3,2,2)
plot(t(1:limit),w_bi(1,1:limit));
title("x-components of the angular velocity in the inertial frame");
xlabel("t(s)");
subplot(3,2,4)
plot(t(1:limit),w_bi(2,1:limit));
title("y-components of the angular velocity in the inertial frame");
xlabel("t(s)");
subplot(3,2,6)
plot(t(1:limit),w_bi(3,1:limit));
title("z-components of the angular velocity in the inertial frame");
xlabel("t(s)");

%% Q2 & Q3 & Bonusses --first run the initial section and run this part--
dt=0.01;
t(1)=1;

limit=1000;

J=zeros(3,3,4);
nu=3.986004418*10^14;
M=[-0.09; 0.001; 0.11];
pos_ecip=zeros(3,25000);
mag_ecip=zeros(3,25000);

J(:,:,1)=[0.018,0,0;0,0.018,0;0,0,0.065]; %Q2
J(:,:,2)=[6.9,0,0;0,7.5,0;0,0,8.4]; %Q3
J(:,:,3)=[6.9,0,0;0,8.4,0;0,0,7.5]; %Bonus-1
J(:,:,4)=[6.9,0.05,0.1;0.05,8.4,0.15;0.15,0.1,7.5]; %Bonus-2


for g=1:4
for s=1:limit
    
    invJ=inv(J(:,:,g));
    t(s)=s*dt;
    q_bi(:,s)=[q_bi(4,s),q_bi(1,s),q_bi(2,s),q_bi(3,s)];%
    A_bi(:,:,s)=quat2dcm(transpose(q_bi(:,s)));
    q_bi(:,s)=[q_bi(2,s),q_bi(3,s),q_bi(4,s),q_bi(1,s)];%

    w_i(:,s)=transpose(A_bi(:,:,s))*w_bi(:,s);
    spinaxis(:,s)=w_i(:,s)/norm(w_i(:,s));

    pos_ecip(1,s)=interp1(1:25000,pos_eci(1,:),t(s),'nearest','extrap');
    pos_ecip(2,s)=interp1(1:25000,pos_eci(2,:),t(s),'nearest','extrap');
    pos_ecip(3,s)=interp1(1:25000,pos_eci(3,:),t(s),'nearest','extrap');

    mag_ecip(1,s)=interp1(1:25000,mag_eci(1,:),t(s),'nearest','extrap');
    mag_ecip(2,s)=interp1(1:25000,mag_eci(2,:),t(s),'nearest','extrap');
    mag_ecip(3,s)=interp1(1:25000,mag_eci(3,:),t(s),'nearest','extrap');


    dist=norm(pos_ecip(:,s));
    n=(A_bi(:,:,s)*-pos_ecip(:,s))/norm(A_bi(:,:,s)*-pos_ecip(:,s));
    N_gg=((3*nu)/dist^3)*cross(n,(J(:,:,g)*n));
    mag_body=A_bi(:,:,s)*(mag_ecip(:,s));
    N_md=cross(M,mag_body)*1e-9;
%     N_gg=[0;0;0];
%     N_md=[0;0;0];
    N_t=N_gg+N_md;

    qk1=0.5*omgof(w_bi(:,s))*q_bi(:,s); %rk4
    wk1=invJ*(N_t-(cross(w_bi(:,s),(J(:,:,g)*w_bi(:,s)))));
    qrk=q_bi(:,s)+qk1*(dt/2);
    wrk=w_bi(:,s)+wk1*(dt/2);
    qk2=0.5*omgof(wrk)*qrk;
    wk2=invJ*(N_t-(cross(wrk,(J(:,:,g)*wrk))));
    qrk=q_bi(:,s)+qk2*(dt/2);
    wrk=w_bi(:,s)+wk2*(dt/2);
    qk3=0.5*omgof(wrk)*qrk;
    wk3=invJ*(N_t-(cross(wrk,(J(:,:,g)*wrk))));
    qrk=q_bi(:,s)+qk3*(dt);
    wrk=w_bi(:,s)+wk3*(dt);
    qk4=0.5*omgof(wrk)*qrk;
    wk4=invJ*(N_t-(cross(wrk,(J(:,:,g)*wrk))));
    q_bi(:,s+1)=q_bi(:,s)+((qk1+2*qk2+2*qk3+qk4)/6)*dt;
    w_bi(:,s+1)=w_bi(:,s)+((wk1+2*wk2+2*wk3+wk4)/6)*dt;

    q_bi(:,s+1)=q_bi(:,s+1)/norm(q_bi(:,s+1));
    

end


figure(g);
subplot(3,2,1)
plot(t(1:limit),spinaxis(1,1:limit));
title("x-components of the spin axis in the inertial frame");
xlabel("t(s)");
ylim([-1 1])
subplot(3,2,3)
plot(t(1:limit),spinaxis(2,1:limit));
title("y-components of the spin axis in the inertial frame");
xlabel("t(s)");
ylim([-1 1])
subplot(3,2,5)
plot(t(1:limit),spinaxis(3,1:limit));
title("z-components of the spin axis in the inertial frame");
xlabel("t(s)");
ylim([-1 1])
hold on;
subplot(3,2,2)
plot(t(1:limit),w_bi(1,1:limit));
title("x-components of the angular velocity in the inertial frame");
xlabel("t(s)");

subplot(3,2,4)
plot(t(1:limit),w_bi(2,1:limit));
title("y-components of the angular velocity in the inertial frame");
xlabel("t(s)");

subplot(3,2,6)
plot(t(1:limit),w_bi(3,1:limit));
title("z-components of the angular velocity in the inertial frame");
xlabel("t(s)");


end


%%

function D = omgof(w)
    D=[-[0,-w(3),w(2);w(3),0,-w(1);-w(2),w(1),0], w(:); -transpose(w(:)), 0];
end