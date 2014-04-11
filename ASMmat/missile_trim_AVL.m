function results = missile_trim_AVL(gb,flightConditions,CLtrim0)

if nargin==0
    clc
    gb = load_configuration(1);
    fc = get_flight_conditions(5000,100);
    results = missile_trim_AVL(gb,fc,0);
    return
end

fc = flightConditions;
cs = control_deflections;
CD0 = gb.CD0;

tol = 1e-4;
iterMax = 5;
err = inf;
iter = 0;
controls = cs.defaults();
x0 = [controls.fin1; controls.alpha];

updCLtrim=0;
if CLtrim0==0
    updCLtrim=1;
end
Cmtrim0 = 0;

while err>=tol && iter<=iterMax
    controls = cs.set_de(controls,x0(2));
    controls.alpha = x0(1);
    aero = missile_analysis_AVL(gb,fc,controls);
    if updCLtrim==1
        CLtrim0 = sqrt(CD0/aero.k);
    end
    CL0 = aero.CL - x0(1)*aero.derivs.CLa - x0(2)*aero.control.CLde;
    Cm0 = aero.Cm - x0(1)*aero.derivs.Cma - x0(2)*aero.control.Cmde;
    CLtrim = CLtrim0 - CL0;
    Cmtrim = Cmtrim0 - Cm0;
    A = [aero.derivs.Cma aero.control.Cmde; aero.derivs.CLa aero.control.CLde];
    b = [Cmtrim; CLtrim];
    xNew = A\b;
    err = norm(x0-xNew);
    x0 = xNew;
    iter = iter +1;
end

results = aero;
end