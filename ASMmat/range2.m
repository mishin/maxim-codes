function [R,aero] = range2(gb, hStart, hEnd, Vstart)

if nargin==0
    clc
    close all
    gb = load_configuration(1);
    tic
    [R,~] = range2(gb,8000, 0, 640/3.6);
    toc
    return
end

CLtrim = 0.9;
Sref = gb.wing.area;
W = gb.weight;
m = gb.Mtot;

fid = fopen('simulation7mod.txt','wt');
fprintf(fid,'time\taltitude\trange\tax\tay\ta\tVx\tVy\tV\tphi\talpha\t');
fprintf(fid,'CL\tCD\telev\tLDmax\n');

lift = .99*W;

dt = 2;

% axes: +X - to the right
% +Y - up
phi = 0;
V  = Vstart;
Vx = Vstart;
Vy = 0;
h = hStart;
R = 0;
t = 0;
i = 0;
figure(1)
hold on
grid on
axis([0,100,0,15])
while h>hEnd
    i = i+1;
    fc = get_flight_conditions(h,V);
    lift = .99*W*cos(phi);
    CLtrim = lift/(fc.dynamicPressure*Sref);
    aero = missile_trim_AVL(gb, fc, CLtrim);
    %lift = aero.CL*fc.dynamicPressure*Sref;
    drag = (gb.CD0+aero.k*aero.CL^2)*fc.dynamicPressure*Sref;
    Fx = lift*sin(phi)-drag*cos(phi);
    Fy = lift*cos(phi)+drag*sin(phi)-W;
    ax = Fx/m;
    ay = Fy/m;
    Vx = Vx + ax*dt;
    Vy = Vy + ay*dt;
    a = sqrt(ax*ax + ay*ay);
    V = sqrt(Vx*Vx + Vy*Vy);
    dh = Vy*dt + ay*dt*dt/2;
    dR = Vx*dt + ax*dt*dt/2;
    dR = dR;
    h = h+dh;
    R = R+dR;
    t = t+dt;
    phi = atan(-Vy/Vx);
    dRange(i) = dR;
    dTime(i)  = dt;
    fprintf(fid,'%.1f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t',t,h,R,ax,ay,a);
    fprintf(fid,'%.4f\t%.4f\t%.4f\t',Vx,Vy,V);
    fprintf(fid,'%.4f\t%.4f\t',rad2deg(phi),aero.alpha);
    fprintf(fid,'%.4f\t%.4f\t%.4f\t%.4f\n',aero.CL,aero.CD,aero.elevator,aero.LD);
    plot(R/1e3,h/1e3,'r*')
end
sum(dRange)
sum(dTime)
fclose(fid);
end