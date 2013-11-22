function results = missile_trim_AVL(gb,flightConditions)

if nargin==0
    clc
    gb = load_configuration(1);
    fc = get_flight_conditions(5000,100);
    results = missile_trim_AVL(gb,fc);
    return
end

fc = flightConditions;
cs = control_deflections;
liftWeightRatio = gb.liftToWeigthRatio;

tol = 1e-4;
iterMax = 5;
err = inf;
iter = 0;
controls = cs.defaults();
x0 = [controls.fin1; controls.alpha];

CLtrim0 = liftWeightRatio*2*gb.mass*9.81/(fc.density*fc.velocity^2*gb.wing.area);
Cmtrim0 = 0;

while err>=tol && iter<=iterMax
    controls = cs.set_de(controls,x0(2));
    controls.alpha = x0(1);
    aero = missile_analysis_AVL(gb,fc,controls);
    CL0 = aero.CL - x0(1)*aero.derivs.CLa - x0(2)*aero.control.CLde;
    Cm0 = aero.Cm - x0(1)*aero.derivs.Cma - x0(2)*aero.control.Cmde;
        CLtrim = CLtrim0 - CL0;
        Cmtrim = Cmtrim0 - Cm0;
    %end
    A = [aero.derivs.Cma aero.control.Cmde; aero.derivs.CLa aero.control.CLde];
    b = [Cmtrim; CLtrim];
    xNew = A\b;
    err = norm(x0-xNew);
    x0 = xNew;
    iter = iter +1;
end

results = aero;
end