function [R,aero] = range(gb, hStart, hEnd, Vstart)

if nargin==0
    clc
    gb = load_configuration(1);
    tic
    [R,~] = range(gb,8500, 1000, 640/3.6);
    toc
    return
end

CLtrim = 0.0;
nSeg = 100;
hStep = -(hEnd-hStart)/nSeg;
Sref = gb.wing.area;
W = gb.weight;
m = gb.Mtot;

dRange = zeros(nSeg,1);
h     = hStart:-hStep:hEnd;
dTime  = zeros(nSeg,1);

phi = 0;
V = Vstart;
Vx = Vstart;
Vy = 0;
% axes: +X - to the right
% +Y - up

for i = 1:nSeg;
%    V   = Vstart+a*dt;
    fc = get_flight_conditions(h(i),V);
    aero = missile_trim_AVL(gb, fc, CLtrim);
    lift = aero.CL*fc.dynamicPressure*Sref;
    drag = (gb.CD0+aero.k*aero.CL^2)*fc.dynamicPressure*Sref;
    Fx = lift*sin(phi)-drag*cos(phi);
    Fy = lift*cos(phi)+drag*sin(phi)-W;
    ax = Fx/m;
    ay = Fy/m;
    a = sqrt(ax*ax + ay*ay);
    dt = (-V+sqrt(V*V+2*ay*hStep))/(2*ay);
    V = V + a*dt;
    Vx = Vx + ax*dt;
    Vy = Vy + ay*dt;
    phi = atan(-Vy/Vx);
    dR = hStep*tan(phi);
    dRange(i) = dR;
    dTime(i) = dt;
    fprintf('h=%.2f\tax=%.2f\tay=%.2f\ta=%.2f\tdt=%.2f\tVx=%.2f\tVy=%.2f\tV=%.2f\tphi=%.2f\talpha=%.2f\tLD=%.2f\n',h(i),ax,ay,a,dt,Vx,Vy,V,rad2deg(phi),aero.alpha,aero.LD)
end

dRange
dTime
plot(cumsum(dRange),h)
end