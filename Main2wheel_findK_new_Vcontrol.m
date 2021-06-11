%%normal motor
%%mw=0.3545;
%%r=0.0335; 
%%mb=29.6-0.3545;
%%l=0.1; 

%%real wheel
mw=3.75*2;
r=0.085;
mb=10;
l=0.6;

iw=0.5*mw*r^2;
ib=mb*l^2;

tw_stall=16;
vmax=24;
tw_0=8;
dseta_0=200*2*pi/60;

k1=tw_stall/vmax;
k2=(k1*vmax-tw_0)/dseta_0;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
c1=(mb+2*(mw+iw/r^2)-mb^2*l^2/(mb*l^2+ib));
c2=-mb^2*l^2*g/(mb*l^2+ib);
c3=2*(1/r+mb*l/(mb*l^2+ib));
c4=(ib+mb*l^2-mb^2*l^2/(mb+2*(mw+iw/r^2)));
c5=(2/r)*(r+mb*l/(mb+2*(mw+iw/r^2)));

A=[0 1 0 0;
   0 -c3*k2/(c1*r) c2/c1 c3*k2/c1;
   0 0 0 1;
   0 c5*k2/(r*c4) mb*l*g/c4 -c5*k2/c4];
B=[0;c3*k1/c1;0;-c5*k1/c4];

C=[1 0 0 0;
   0 0 1 0];
D=[0 ;0];
Q=[1000 0 0 0;
   0 1 0 0;
   0 0 100 0;
   0 0 0 1];
R=0.1;
K = lqr(A,B,Q,R)

Ac = [(A-B*K)];
Bc = [B];
Cc = [C];
Dc = [D];

states = {'x' 'x_dot' 'seta' 'seta_dot'};
inputs = {'tw'};
outputs = {'x' 'seta'};

sys_cl = ss(Ac,Bc,Cc,Dc,'statename',states,'inputname',inputs,'outputname',outputs);

t = 0:0.01:10;
r =0*ones(size(t));
r(1,1)=1;
[y,t,x]=lsim(sys_cl,r,t);
[AX,H1,H2] = plotyy(t,y(:,1),t,y(:,2),'plot');
set(get(AX(1),'Ylabel'),'String','wheel position (m)')
set(get(AX(2),'Ylabel'),'String','pendulum angle (radians)')
title('Step Response with LQR Control')


