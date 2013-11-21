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

tol = 1e-6;
iterMax = 10;
err = inf;
iter = 0;
controls = cs.defaults();
x0 = [controls.fin1; controls.alpha];
aero = missile_analysis_AVL(gb,fc,controls);
CLtrim0 = liftWeightRatio*2*gb.mass*9.81/(fc.density*fc.velocity^2*gb.wing.area) - aero.CL;
Cmtrim0 = - aero.Cm;

while err>=tol && iter<=iterMax
    if iter>0
        aero = missile_analysis_AVL(gb,fc,controls);
    end
    A = [aero.control.Cmde aero.derivs.Cma; aero.control.CLde aero.derivs.CLa];
    b = [Cmtrim0; CLtrim0];
    xnew = A\b;
    err = norm(x0-xnew);
    x0 = xnew;
    iter = iter +1;
    controls = cs.set_de(controls,x0(1));
    controls.alpha = x0(2);
    %CL(iter) = aero.CL;
    %plot(1:iter,CL,'bo-');
end

results = aero;
end