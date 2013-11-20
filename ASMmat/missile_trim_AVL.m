function results = missile_trim_AVL(gb,flightConditions)

if nargin==0
    gb = load_configuration();
    fc = get_flight_conditions(5000,100);
    results = missile_trim_AVL(gb,fc);
    return
end

fc = flightConditions;

errCm = 1e-6;
iterMax = 5;
err = inf;
iter = 0;

while err>=errCm || iter<=iterMax
    aero = missile_analysis_AVL(gb,fc);
end

end