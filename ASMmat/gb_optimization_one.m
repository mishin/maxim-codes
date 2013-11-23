function xOptimum = gb_optimization_one()
clc
clear all
%% input
x0 = [1.3995 0.0889 1.0 0.0 0.85 0.12 0.0762 1.0 1.6];
xL = [1.2000 0.0600 0.6 0.0 0.6 0.07 0.0600 0.6 1.4];
xU = [1.7000 0.1200 1.0 30. 1.2 0.15 0.0900 1.0 1.7];
hStart = 8500;
hEnd = 1000;
Vstart = 180;

% constraints
SMmin = 0.05;
SMmax = 0.30;
ClbetaMax = -0.02; %-0.03
CnbetaMin = 0.03; %0.08
CnbetaMax = 0.28;
CmalphaMax = 0;
rangeMin = 70;
alphaTrimMin = -10;
alphaTrimMax = 10;
elevatorTrimMin = -15;
elevatorTrimMax = 15;


%% calculation
    function [f,g] = run_analysis(x)
        gb = load_configuration(1);
        gb.wing.span = x(1);
        gb.wing.secChords = [x(2), x(2)*x(3)];
        gb.wing.sweep = x(4);
        gb.wing.location(1) = x(5);
        gb.fin.halfSpan = x(6);
        gb.fin.secChords = [x(7), x(7)*x(8)];
        gb.fin.locationX = x(9);
        gb = geometry_analysis(gb);
        gb.massData = weight_analysis(gb);
        gb.CG = [gb.massData.Xcg 0 0];
        [Range, aero] = range(gb,hStart,hEnd,Vstart);
        f = -Range;
        g(1) = SMmin - aero.SM;
        g(2) = aero.SM - SMmax;
        %g(3) = aero.derivs.Clb - ClbetaMax;
        g(3) = -1;
        g(4) = CnbetaMin - aero.derivs.Cnb;
        g(5) = aero.derivs.Cnb - CnbetaMax;
        g(6) = aero.derivs.Cma - CmalphaMax;
        g(7) = rangeMin - Range;
        g(8) = alphaTrimMin - aero.alpha;
        g(9) = aero.alpha - alphaTrimMax;
        g(10) = elevatorTrimMin - aero.elevator;
        g(11) = aero.elevator - elevatorTrimMax;
        fprintf('%.4f  %.4f  %.4f  %.4f  %.4f  %.4f\n' ,aero.SM, aero.derivs.Clb, aero.derivs.Cnb, aero.derivs.Cma,aero.alpha,aero.elevator);
    end

    function f = objective(x)
        [obj,~] = run_analysis(x);
        f = obj;
    end
    function [g,h] = constraints(x)
        [~,cnstr] = run_analysis(x);
        g = cnstr;
        h = [];
    end



options = optimset('GradObj','off','Algorithm','active-set','FinDiffType','forward','GradConstr','off','TolX',1e-6,'MaxIter',30,'Display','iter');
xOptimum = fmincon(@objective,x0,[],[],[],[],xL,xU,@constraints,options);

R0 = objective(x0);
R1 = objective(xOptimum);
constraints(xOptimum)
fprintf('Range baseline: %.4f km\nRange optimum%.4f km\n',R0, R1);

end