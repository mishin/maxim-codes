function fc = get_flight_conditions(altitude, velocity)

if nargin<2
    fc = get_flight_conditions(0,100);
    return
end

%% datcom required input
fc.alpha = [-10,10];
fc.beta = [-10,10];
fc.deflection = [0, 0, 0, 0];

%% standard atmosphere
fc.altitude = altitude;
P0 = 101325;
T0 = 288.15;
G0 = 9.80665;
R = 287.04;
gamma = 1.4;
dT = 0;

if altitude<11000
    T = T0 - 6.5 * altitude / 1000.0 + dT;
    P = P0 * (1 - 0.0065 * altitude / T0)^5.2561;
else
    T11 = T0 - 6.5 * 11000 / 1000;
    P11 = P0 * (1 - 0.0065 * 11000 / T0)^5.2561;
    T = T11 + dT;
    P = P11 * exp(-G0/(R*T11) * (altitude - 11000));
end

fc.density = P/(R*T);
fc.temperature = T;
fc.pressure = P;
fc.soundSpeed = (gamma*R*T)^0.5;

if velocity<1
    fc.Mach = velocity;
    fc.velocity = fc.Mach*fc.soundSpeed;
else
    fc.velocity = velocity;
    fc.Mach = velocity/fc.soundSpeed;
end

fc.dynamicPressure = fc.density*fc.velocity^2/2;

end