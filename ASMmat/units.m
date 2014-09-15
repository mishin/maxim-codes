function answer = units( x, unit1, unit2 );
%
% y = units(x,'unit1','unit2')  - Converts x from physical unit1 to unit2
% example: units(1,'in','mm') = 25.4
%
% The following units are supported:
% Acceleration: m/s^2, cm/s^2, mm/s^2, ft/s^2, in/s^2, G
% Angle: rad, deg, rev
% Area: km^2, m^2, cm^2, mm^2, ym^2 (square-micrometer), sq-mile, sq-yd, sq-ft, sq-in,
%    acres, ha, ar
% Area Moment of Inertia: m^4, cm^4, mm^4, ft^4, in^4
% Density: t/m^3 (metric), kg/m^3, g/cm^3, g/mm^3, lbs/ft^3 lbs/in^3, lbs/galUS, lbs/galUK
% Energy, Work, Torque: GJ, MJ, kJ, J, mJ, Nm, Ncm, Nmm, kWh, Wh, Ws, lb-ft, lb-in, oz-in,Btu,
%    Btu, cal, kcal, eV
% Force: MN, kN, N, dyne, lbf, kip
% Fuel Consumption: l/100km, miles/galUS
% Frequency, Angular Velocity: GHz, MHz, kHz, Hz, 1/s, 1/min, 1/h, rad/s, deg/s, rpm
% Length: km, m, dm, cm, mm, ym (micrometer), nm, mile, yard, ft, in, mill,
%    Angstrom, light-year, parsec
% Mass: t (metric), tUS, tUK, kg, g, mg, yg (microgram), ng, lbs, oz
% Mass Moment of Inertia: kg*m^2, kg*cm^2, kg*mm^2, g*m^2, g*cm^2, g*mm^2, lb*ft^2, lb*in^2
% Power: GW, MW, kW, W, mW, hp, Btu/h, Btu/s, kcal/h, J/h
% Pressure, Stress: GPa, MPa, kPa, hPa, Pa, bar, mbar, atm, dyne/cm^2, ksi, psi, mmHg, mmH2O
% Strain: m/m, mm/m, ym/m (micrometer/m), nm/m, %, o/oo, in/in, mill/in
% Stress Intensity Factor: MPa*m^1/2, MPa*mm^1/2, ksi*in^1/2, psi*in^1/2
% Temperature: degK, degC, degF, degR
% Time: yr (365 days), mth (30 days), wk, day, hr, min, s, ms, ys (microsecond), ns
% Velocity: km/s, km/h, m/s, cm/s, mm/s, m/min, mm/min, mps, mph, ft/s, in/s, ft/min,
%    in/min, Mach, knots
% Viscosity: Ns/m^2, poise, centipoise, lbfs/sq-ft
% Volume: km^3, m^3, cm^3, mm^3, ym^3, cu-mile, cu-ft, cu-in, l, cl, ml, galUS, galUK,
%    pint (liquid US), quart (liquid US), fl-oz (liquid US)
% Metric Prefixes: Yotta, Zetta, Exa, Peta, Terra, Giga, Mega, Myria, kilo, hecto, 1
%     deci, centi, milli, micro, nano, pico, femto, atto, zepto, yocto
% -------------------------------------------------------------------------
% Henning Ressing, PhD
% The University of British Columbia, Mechanical Engineering
% Updated and modified by W. Scott Richardson
% B.S.E. Aerospace Engineering, University of Michigan - Ann Arbor
% Last updated: July 27, 2004



err = 0;
unitpick=0;
if unitpick==0
   switch unit1
      % Length
   case 'km'
      unitpick=1;
      x = x.*1000; 
   case 'm'
      unitpick=1;
      x = x;
   case 'dm'
      unitpick=1;
      x = x./10;
   case 'cm'
      unitpick=1;
      x = x./100;
   case 'mm'
      unitpick=1;
      x = x./1000;
   case 'ym'
      unitpick=1;
      x = x./1e6;
   case 'nm'
      unitpick=1;
      x = x./1e9;
   case 'mile'
      unitpick=1;
      x = x.*1609.344;
   case 'yard'
      unitpick=1;
      x = x.*0.9144;
   case 'ft'
      unitpick=1;
      x = x.*0.3048;  
   case 'in'
      unitpick=1;
      x = x.*0.0254;
   case 'mill'
      unitpick=1;
      x = x.*(25.4e-6);
   case 'Angstrom'
      unitpick=1;
      x = x./1e10;
   case 'light-year'
      unitpick=1;
      x = x.*9.46075309081898e15;
   case 'parsec'
      unitpick=1;
      x = x.*3.262*9.46075309081898e15;
   end
end

if unitpick==0
   switch unit1
      % Time
   case 'yr'
      unitpick=2;
      x = x.*(60*60*24*365);  
   case 'mth'
      unitpick=2;
      x = x.*(60*60*24*30);
   case 'wk'
      unitpick=2;
      x = x.*(60*60*24*7);
   case 'day'
      unitpick=2;
      x = x.*(60*60*24);
   case 'hr'
      unitpick=2;
      x = x.*(60*60);
   case 'min'
      unitpick=2;
      x = x.*60;
   case 's'
      unitpick=2;
      x = x;
   case 'ms';
      unitpick=2;
      x = x./1000;
   case 'ys';
      unitpick=2;
      x = x./1e6;
   case 'ns';
      unitpick=2;
      x = x./1e9;
   end
end

if unitpick==0
   switch unit1
      % Mass
   case 't'
      unitpick=3;
      x = x.*1000;
   case 'tUS'
      unitpick=3;
      x = x.*907.18474;
   case 'tUK'
      unitpick=3;
      x = x.*1016.0469088;
   case 'kg'
      unitpick=3;
      x = x;  
   case 'g'
      unitpick=3;
      x = x./1000;
   case 'mg'
      unitpick=3;
      x = x./1e6;
   case 'yg'
      unitpick=3;
      x = x./1e9;
   case 'ng'
      unitpick=3;
      x = x./1e12;
   case 'lbs'
      unitpick=3;
      x = x.*0.4535924;
   case 'oz'
      unitpick=3;
      x = x.*0.03110348;
   end
end

if unitpick==0
   switch unit1
      % Force
   case 'MN'
      unitpick=4;
      x = x.*1e6;  
   case 'kN'
      unitpick=4;
      x = x.*1000;  
   case 'N'
      unitpick=4;
      x = x;
   case 'dyne'
      unitpick=4;
      x = x./1e5;
   case 'lbf'
      unitpick=4;
      x = x.*4.448222;
   case 'kip'
      unitpick=4;
      x = x.*4448.222
   end
end

if unitpick==0
   switch unit1
      % Pressure, Stress
   case 'GPa'
      unitpick=5;
      x = x.*1e9;  
   case 'MPa'
      unitpick=5;
      x = x.*1e6;
   case 'N/mm^2'
      unitpick=5;
      x = x.*1e6;
   case 'kPa'
      unitpick=5;
      x = x.*1000;
   case 'hPa'
      unitpick=5;
      x = x.*100;
   case 'Pa'
      unitpick=5;
      x = x;
   case 'bar'
      unitpick=5;
      x = x.*1e5;
   case 'mbar'
      unitpick=5;
      x = x.*100;
   case 'dyne/cm^2'
      unitpick=5;
      x = x./10;
   case 'atm'
      unitpick=5;
      x = x.*1.01325e5;
   case 'ksi'
      unitpick=5;
      x = x.*6894757;
   case 'psi'
      unitpick=5;
      x = x.*6894.757;
   case 'mmHg'
      unitpick=5;
      x = x.*133.322;
   case 'mmH2O'
      unitpick=5;
      x = x.*9.80665;
   end
end

if unitpick==0
   switch unit1
      % Temperature
   case 'degK'
      unitpick=6;
      x = x;
   case 'degC'
      unitpick=6;
      x = x + 273.15;
   case 'degF'
      unitpick=6;
      x = 5/9.*(x-32) + 273.15;
   case 'degR'
      unitpick=6;
      x = 5/9.*(x-491.67) + 273.15;
   end
end

if unitpick==0
   switch unit1
      % Work, Energy, Torque
   case 'GJ'
      unitpick=7;
      x = x.*10^9;
   case 'MJ'
      unitpick=7;
      x = x.*10^6;
   case 'kJ'
      unitpick=7;
      x = x.*1000;
   case 'J'
      unitpick=7;
      x = x;
   case 'mJ'
      unitpick=7;
      x = x./1000;
   case 'Nm'
      unitpick=7;
      x = x;
   case 'Ncm'
      unitpick=7;
      x = x./100;
   case 'Nmm'
      unitpick=7;
      x = x./1000;
   case 'lb-ft'
      unitpick=7;
      x = x.*1.355818;
   case 'lb-in'
      unitpick=7;
      x = x.*0.112984833;
   case 'oz-in'
      unitpick=7;
      x = x.*7.061552e-3;
   case 'Btu'
      unitpick=7;
      x = x.*1055.056;
   case 'kWh'
      unitpick=7;
      x = x.*3.6e6;
   case 'Wh'
      unitpick=7;
      x = x.*3600;
   case 'Ws'
      unitpick=7;
      x = x;
   case 'cal'
      unitpick=7;
      x = x.*4.1868;
   case 'kcal'
      unitpick=7;
      x = x.*4186.8;
   case 'eV'
      unitpick=7;
      x = x.*1.60218e-19;
   end
end

if unitpick==0
   switch unit1
      % Power
   case 'GW'
      unitpick=8;
      x = x.*1e9;
   case 'MW'
      unitpick=8;
      x = x.*1e6;
   case 'kW'
      unitpick=8;
      x = x.*1e3;
   case 'W'
      unitpick=8;
      x = x;
   case 'mW'
      unitpick=8;
      x = x./1000;
   case 'hp'
      unitpick=8;
      x = x.*745.6999;
   case 'Btu/h'
      unitpick=8;
      x = x.*0.2930711;
   case 'Btu/s'
      unitpick=8;
      x = x.*1055.056;
   case 'kcal/h'
      unitpick=8;
      x = x.*1.163;
   case 'J/h'
      unitpick=8;
      x = x./3600;
   end
end

if unitpick==0
   switch unit1
      % Velocity
   case 'km/h'
      unitpick=9;
      x = x.*(1000/3600);
   case 'km/s'
      unitpick=9;
      x = x.*1000;
   case 'm/s'
      unitpick=9;
      x = x;
   case 'cm/s'
      unitpick=9;
      x = x./100;
   case 'mm/s'
      unitpick=9;
      x = x./1000;
   case 'm/min'
      unitpick=9;
      x = x./60;
   case 'mm/min'
      unitpick=9;
      x = x./60000;
   case 'mph'
      unitpick=9;
      x = x.*(1609.344/3600);
   case 'mps'
      unitpick=9;
      x = x.*1609.344;
   case 'ft/s'
      unitpick=9;
      x = x.*0.3048;
   case 'in/s'
      unitpick=9;
      x = x.*0.0254;
   case 'ft/min'
      unitpick=9;
      x = x.*(0.3048/60);
   case 'in/min'
      unitpick=9;
      x = x.*(0.0254/60);
   case 'Mach'
      unitpick=9;
      x = x.*331.5;
   case 'knots'
      unitpick=9;
      x = x.*0.514444;
   end
end

if unitpick==0
   switch unit1
      % Acceleration
   case 'm/s^2'
      unitpick=10;
      x = x;
   case 'cm/s^2'
      unitpick=10;
      x = x./100;
   case 'mm/s^2'
      unitpick=10;
      x = x./1000;
   case 'ft/s^2'
      unitpick=10;
      x = x.*0.3048;
   case 'in/s^2'
      unitpick=10;
      x = x.*0.0254;
   case 'G'
      unitpick=10;
      x = x.*9.81;
   end
end

if unitpick==0
   switch unit1
      % Area
   case 'km^2'
      unitpick=11;
      x = x.*1e6;
   case 'm^2'
      unitpick=11;
      x = x;
   case 'cm^2'
      unitpick=11;
      x = x./1e4;
   case 'mm^2'
      unitpick=11;
      x = x./1e6;
   case 'ym^2'
      unitpick=11;
      x = x./1e12;
   case 'sq-mile'
      unitpick=11;
      x = x.*1609.344^2;
   case 'sq-yd'
      unitpick=11;
      x = x.*0.83612736;
   case 'sq-ft'
      unitpick=11;
      x = x.*0.3048^2;
   case 'sq-in'
      unitpick=11;
      x = x.*0.0254^2;
   case 'acres'
      unitpick=11;
      x = x.*4046.8564224;
   case 'ha'
      unitpick=11;
      x = x.*1e4;
   case 'ar'
      unitpick=11;
      x = x.*100;
   end
end

if unitpick==0
   switch unit1
      % Volume
   case 'km^3'
      unitpick=12;
      x = x.*10^9;
   case 'm^3'
      unitpick=12;
      x = x;
   case 'cm^3'
      unitpick=12;
      x = x./10^6;
   case 'mm^3'
      unitpick=12;
      x = x./10^9;
   case 'ym^3'
      unitpick=12;
      x = x./10^18;
   case 'cu-mile'
      unitpick=12;
      x = x.*1609.344^3;
   case 'cu-yd'
      unitpick=12;
      x = x.*0.764554857984;
   case 'cu-ft'
      unitpick=12;
      x = x.*0.3048^3;
   case 'cu-in'
      unitpick=12;
      x = x.*0.0254^3;
   case 'l'
      unitpick=12;
      x = x./1000;
   case 'cl'
      unitpick=12;
      x = x./1e6;
   case 'ml'
      unitpick=12;
      x = x./10^9;
   case 'galUS'
      unitpick=12;
      x = x.*3.785412e-3;
   case 'galUK'
      unitpick=12;
      x = x.*4.54609e-3; 
   case 'pint'
      unitpick=12;
      x = x.*4.73176473e-4;
   case 'quart'
      unitpick=12;
      x = x.*9.46352946e-4;
   case 'fl-oz'
      unitpick=12;
      x = x.*2.95735295625e-5;
   end
end

if unitpick==0
   switch unit1
      % Density
   case 't//m^3'
      unitpick=13;
      x = x.*1000;
   case 'kg/m^3'
      unitpick=13;
      x = x;
   case 'g/cm^3'
      unitpick=13;
      x = x.*1000;
   case 'g/mm^3'
      unitpick=13;
      x = x.*1e6;
   case 'kg/l'
      unitpick=13;
      x = x.*1000;
   case 'lbs/ft^3'
      unitpick=13;
      x = x.*(0.4535924/0.3048^3);
   case 'lbs/in^3'
      unitpick=13;
      x = x.*(0.4535924/0.0254^3);
   case 'lbs/galUS'
      unitpick=13;
      x = x.*119.826427;
   case 'lbs/galUK'
      unitpick=13;
      x = x.*99.776373;
   end
end

if unitpick==0
   switch unit1
      % Mass Moment of Inertia
   case 'kg*m^2'
      unitpick=14;
      x = x;
   case 'kg*cm^2'
      unitpick=14;
      x = x./1e4;
   case 'kg*mm^2'
      unitpick=14;
      x = x./1e6;
   case 'g*m^2'
      unitpick=14;
      x = x./1000;
   case 'g*cm^2'
      unitpick=14;
      x = x./1e7;
   case 'g*mm^2'
      unitpick=14;
      x = x./1e9;
   case 'lb*ft^2'
      unitpick=14;
      x = x.*0.04214011;
   case 'lb*in^2'
      unitpick=14;
      x = x.*0.2926397e-3;
   end
end

if unitpick==0
   switch unit1
      % Area Moment of Inertia
   case 'm^4'
      unitpick=15;
      x = x;
   case 'cm^4'
      unitpick=15;
      x = x./1e8;
   case 'mm^4'
      unitpick=15;
      x = x./1e12;
   case 'ft^4'
      unitpick=15;
      x = x./(0.3048^4);
   case 'in^4'
      unitpick=15;
      x = x./(0.0254^4);
   end
end

if unitpick==0
   switch unit1
      % Frequency / Angular Velocity
   case 'GHz'
      unitpick=16;
      x = x.*1e9
   case 'MHz'
      unitpick=16;
      x = x.*1e6
   case 'kHz'
      unitpick=16;
      x = x.*1e3
   case 'Hz'
      unitpick=16;
      x = x;
   case '1/min'
      unitpick=16;
      x = x./60;
   case '1/h'
      unitpick=16;
      x = x./3600;
   case 'rad/s'
      unitpick=16;
      x = x./(2*pi);
   case 'deg/s'
      unitpick=16;
      x = x./360;
   case 'rpm'
      unitpick=16;
      x = x./60;
   end
end

if unitpick==0
   switch unit1
      % Angle
   case 'rad'
      unitpick=17;
      x = x;
   case 'deg'
      unitpick=17;
      x = x.*(pi/180);
   case 'rev'
      unitpick=17;
      x = x.*(2*pi);
   end
end

if unitpick==0
   switch unit1
      % Stress Intensity Factor
   case 'MPa*m^1/2'
      unitpick=18;
      x = x;
   case 'MPa*mm^1/2'
      unitpick=18;
      x = x./sqrt(1000);
   case 'ksi*in^1/2'
      unitpick=18;
      x = x.*6.894757.*sqrt(0.0254);
   case 'psi*in^1/2'
      unitpick=18;
      x = x.*6894.757.*sqrt(0.0254);
   end
end

if unitpick==0
   switch unit1
      % Fuel Consumption
   case 'l/100km'
      unitpick=19;
      x = x;
   case 'miles/galUS'
      unitpick=19;
      x = 235.214596754951./x;
   end
end

if unitpick==0
   switch unit1
      % Viscosity
   case 'Ns/m^2'
      unitpick=20;
      x = x;
   case 'poise'
      unitpick=20;
      x = x./10;
   case 'centipoise'
      unitpick=20;
      x = x./1000;
   case 'lbfs/sq-ft'
      unitpick=20;
      x = x./0.02089;
   end
end

if unitpick==0
   switch unit1
      % Strain
   case 'm/m'
      unitpick=21;
      x = x;
   case 'mm/m'
      unitpick=21;
      x = x./1e3;
   case 'ym/m'
      unitpick=21;
      x = x./1e6;
   case 'nm/m'
      unitpick=21;
      x = x./1e9;
   case '%'
      unitpick=21;
      x = x./1e2;
   case 'o/oo'
      unitpick=21;
      x = x./1e3;
   case 'in/in'
      unitpick=21;
      x = x;
   case 'mill/in'
      unitpick=21;
      x = x./1e3';
   end
end

if unitpick==0
   switch unit1
      % Metric Prefixes
   case 'Yotta'
      unitpick=22;
      x = x.*1e24;
   case 'Zetta'
      unitpick=22;
      x = x.*1e21;
   case 'Exa'
      unitpick=22;
      x = x.*1e18;
   case 'Peta'
      unitpick=22;
      x = x.*1e15;
   case 'Tera'
      unitpick=22;
      x = x.*1e12;
   case 'Giga'
      unitpick=22;
      x = x.*1e9;
   case 'Mega'
      unitpick=22;
      x = x.*1e6;
   case 'Myria'
      unitpick=22;
      x = x.*1e5;
   case 'kilo'
      unitpick=22;
      x = x.*1e3;
   case 'hecto'
      unitpick=22;
      x = x.*1e2;
   case '1'
      unitpick=22;
      x = x;
   case 'deci'
      unitpick=22;
      x = x.*1e-1;
   case 'centi'
      unitpick=22;
      x = x.*1e-2;
   case 'milli'
      unitpick=22;
      x = x.*1e-3;
   case 'micro'
      unitpick=22;
      x = x.*1e-6;
   case 'nano'
      unitpick=22;
      x = x.*1e-9;
   case 'pico'
      unitpick=22;
      x = x.*1e-12;
   case 'femto'
      unitpick=22;
      x = x.*1e-15;
   case 'atto'
      unitpick=22;
      x = x.*1e-18;
   case 'zepto'
      unitpick=22;
      x = x.*1e-21;
   case 'yocto'
      unitpick=22;
      x = x.*1e-24;
   otherwise
      disp('Unsupported unit - please restart program!');
      err=1;
   end
end

% ----------------------------------------------------------

if ~err
   if unitpick==1
      switch unit2
         % Length
      case 'km'
         x = x./1000;  
      case 'm'
         x = x;
      case 'dm'
         x = x.*10;
      case 'cm'
         x = x.*100;
      case 'mm'
         x = x.*1000;
      case 'ym'
         x = x.*1e6;
      case 'nm'
         x = x.*1e9;
      case 'mile'
         x = x./1609.344;
      case 'yard'
         x = x./0.9144;
      case 'ft'
         x = x./0.3048;  
      case 'in'
         x = x./0.0254;
      case 'mill'
         x = x./25.4e-6;
      case 'Angstrom'
         x = x.*1e10;
      case 'light-year'
         x = x./9.46075309081898e15;
      case 'parsec'
         x = x.*3.262*9.46075309081898e15;
      otherwise
         disp('Inconsistent unit - please check units and try again!')
         err=2;
      end
   end
end
if ~err
   if unitpick==2
      switch unit2
         % Time
      case 'yr'
         x = x./(60*60*24*365);  
      case 'mth'
         x = x./(60*60*24*30);
      case 'wk'
         x = x./(60*60*24*7);
      case 'day'
         x = x./(60*60*24);
      case 'hr'
         x = x./(60*60);
      case 'min'
         x = x./60;
      case 's'
         x = x;
      case 'ms';
         x = x.*1000;
      case 'ys';
         x = x.*1e6;
      case 'ns';
         x = x.*1e9;
      otherwise
         disp('Inconsistent unit - please check units and try again!')
         err=2;
      end
   end
end
if ~err
   if unitpick==3
      switch unit2
         % Mass
      case 't'
         x = x./1000;
      case 'tUS'
         x = x./907.18474;
      case 'tUK'
         x = x./1016.0469088;
      case 'kg'
         x = x;  
      case 'g'
         x = x.*1000;
      case 'mg'
         x = x.*1e6;
      case 'yg'
         x = x.*1e9;
      case 'ng'
         x = x.*1e12;
      case 'lbs'
         x = x./0.4535924;
      case 'oz'
         x = x./0.03110348;
      otherwise
         disp('Inconsistent unit - please check units and try again!')
         err=2;
      end
   end
end
if ~err
   if unitpick==4
      switch unit2
         % Force
      case 'MN'
         x = x./1e6;  
      case 'kN'
         x = x./1000;  
      case 'N'
         x = x;
      case 'dyne'
         x = x.*1e5;
      case 'lbf'
         x = x./4.44822;
      case 'kip'
         x = x./4448.222
      otherwise
         disp('Inconsistent unit - please check units and try again!')
         err=2;
      end
   end
end
if ~err
   if unitpick==5
      switch unit2
         % Pressure, Stress
      case 'GPa'
         x = x./10^9;  
      case 'MPa'
         x = x./10^6;  
      case 'N/mm^2'
         x = x./10^6;  
      case 'kPa'
         x = x./1000;
      case 'hPa'
         x = x./100;  
      case 'Pa'
         x = x;
      case 'bar'
         x = x./1e5;
      case 'mbar'
         x = x./100;
      case 'dyne/cm^2'
         x = x.*10;
      case 'atm'
         x = x./1.01325e5;
      case 'ksi'
         x = x./6894757;
      case 'psi'
         x = x./6894.757;
      case 'mmHg'
         x = x./133.322;
      case 'mmH2O'
         x = x./9.80665;
      otherwise
         disp('Inconsistent unit - please check units and try again!')
         err=2;
      end
   end
end
if ~err
   if unitpick==6
      switch unit2
         % Temperature
      case 'degK'
         x = x;
      case 'degC'
         temp = size(c);
         x = c - 273.15;
      case 'degF'
         temp = size(c);
         x = 9/5.*(c-273.15) + 32;
      case 'degR'
         temp = size(c);
         x = 9/5.*(c-273.15) + 491.67;
      otherwise
         disp('Inconsistent unit - please check units and try again!')
         err=2;
      end
   end
end
if ~err
   if unitpick==7
      switch unit2
         % Work, Energy, Torque
      case 'GJ'
         x = x./10^9;
      case 'MJ'
         x = x./10^6;
      case 'kJ'
         x = x./1000;
      case 'J'
         x = x;
      case 'mJ'
         x = x.*1000;
      case 'Nm'
         x = x;
      case 'Ncm'
         x = x.*100;
      case 'Nmm'
         x = x.*1000;
      case 'lb-ft'
         x = x./1.355818;
      case 'lb-in'
         x = x./0.112984833;
      case 'oz-in'
         x = x./0.007061552;
      case 'Btu'
         x = x./1055.056;
      case 'kWh'
         x = x./3.6e6;
      case 'Wh'
         x = x./3600;
      case 'Ws'
         x = x;
      case 'cal'
         x = x./4.1868;
      case 'kcal'
         x = x./4186.8;
      case 'eV'
         x = x./1.60218e-19;
      otherwise
         disp('Inconsistent unit - please check units and try again!')
         err=2;
      end
   end
end
if ~err
   if unitpick==8
      switch unit2
         % Power
      case 'GW'
         x = x./1e9;
      case 'MW'
         x = x./1e6;
      case 'kW'
         x = x./1000;
      case 'W'
         x = x;
      case 'mW'
         x = x.*1000;
      case 'hp'
         x = x./745.6999;
      case 'Btu/h'
         x = x./0.2930711;
      case 'Btu/s'
         x = x./1055.056;
      case 'kcal/h'
         x = x./1.163;
      case 'J/h'
         x = x.*3600;
      otherwise
         disp('Inconsistent unit - please check units and try again!')
         err=2;
      end
   end
end
if ~err
   if unitpick==9
      switch unit2
         % Velocity
      case 'km/h'
         x = x./(1000/3600);
      case 'km/s'
         x = x./1000;
      case 'm/s'
         x = x;
      case 'cm/s'
         x = x.*100;
      case 'mm/s'
         x = x.*1000;
      case 'm/min'
         x = x.*60;
      case 'mm/min'
         x = x.*60000;
      case 'mph'
         x = x./(1609.344/3600);
      case 'mps'
         x = x./1609.344;
      case 'ft/s'
         x = x./0.3048;
      case 'in/s'
         x = x./0.0254;
      case 'ft/min'
         x = x./(0.3048*60);
      case 'in/min'
         x = x./(0.0254*60);
      case 'Mach'
         x = x./331.5;
      case 'knots'
         x = x./0.514444;
      otherwise
         disp('Inconsistent unit - please check units and try again!')
         err=2;
      end
   end
end
if ~err
   if unitpick==10
      switch unit2
         % Acceleration
      case 'm/s^2'
         x = x;
      case 'cm/s^2'
         x = x.*100;
      case 'mm/s^2'
         x = x.*1000;
      case 'ft/s^2'
         x = x./0.3048;
      case 'in/s^2'
         x = x./0.0254;
      case 'G'
         x = x./9.81;
      otherwise
         disp('Inconsistent unit - please check units and try again!')
         err=2;
      end
   end
end
if ~err
   if unitpick==11
      switch unit2
         % Area
      case 'km^2'
         x = x./1e6;
      case 'm^2'
         x = x;
      case 'cm^2'
         x = x.*10^4;
      case 'mm^2'
         x = x.*10^6;
      case 'ym^2'
         x = x.*10^12;
      case 'sq-mile'
         x = x./1609.344^2;
      case 'sq-yd'
         x = x./0.83612736;
      case 'sq-ft'
         x = x./0.3048^2;
      case 'sq-in'
         x = x./0.0254^2;
      case 'acres'
         x = x./4046.8564224;
      case 'ha'
         x = x./1e4;
      case 'ar'
         x = x./100;
      otherwise
         disp('Inconsistent unit - please check units and try again!')
         err=2;
      end
   end
end
if ~err
   if unitpick==12
      switch unit2
         % Volume
      case 'km^3'
         x = x./10^9;
      case 'm^3'
         x = x;
      case 'cm^3'
         x = x.*10^6;
      case 'mm^3'
         x = x.*10^9;
      case 'ym^3'
         x = x.*10^18;
      case 'cu-mile'
         x = x./1609.344^3;
      case 'cu-ft'
         x = x./0.3048^3;
      case 'cu-yd'
         x = x./0.764554857984;
      case 'cu-in'
         x = x./0.0254^3;
      case 'l'
         x = x.*1000;
      case 'cl'
         x = x.*1e6;
      case 'ml'
         x = x.*10^9;
      case 'galUS'
         x = x./3.785412e-3;
      case 'galUK'
         x = x./4.54609e-3; 
      case 'pint'
         x = x./4.73176473e-4;
      case 'quart'
         x = x./9.46352946e-4;
      case 'fl-oz'
         x = x./2.95735295625e-5;
      otherwise
         disp('Inconsistent unit - please check units and try again!')
         err=2;
      end
   end
end
if ~err
   if unitpick==13
      switch unit2
         % Density
      case 't/m^3'
         x = x./1000;
      case 'kg/m^3'
         x = x;
      case 'g/cm^3'
         x = x./1000;
      case 'g/mm^3'
         x = x./1e6;
      case 'kg/l'
         x = x./1000;
      case 'lbs/ft^3'
         x = x./(0.4535924/0.3048^3);
      case 'lbs/in^3'
         x = x./(0.4535924/0.0254^3);
      case 'lbs/galUS'
         x = x./119.826427;
      case 'lbs/galUK'
         x = x./99.776373;
      otherwise
         disp('Inconsistent unit - please check units and try again!')
         err=2;
      end
   end
end
if ~err
   if unitpick==14
      switch unit2
         % Mass Moment of Inertia
      case 'kg*m^2'
         x = x;
      case 'kg*cm^2'
         x = x.*1e4;
      case 'kg*mm^2'
         x = x.*1e6;
      case 'g*m^2'
         x = x.*1000;
      case 'g*cm^2'
         x = x.*1e7;
      case 'g*mm^2'
         x = x.*1e9;
      case 'lb*ft^2'
         x = x./0.04214011;
      case 'lb*in^2'
         x = x./0.2926397e-3;
      otherwise
         disp('Inconsistent unit - please check units and try again!')
         err=2;
      end
   end
end
if ~err
   if unitpick==15
      switch unit2
         % Areal Moment of Inertia
      case 'm^4'
         x = x;
      case 'cm^4'
         x = x.*1e8;
      case 'mm^4'
         x = x.*1e12;
      case 'ft^4'
         x = x.*(0.3048^4);
      case 'in^4'
         x = x.*(0.0254^4);
      otherwise
         disp('Inconsistent unit - please check units and try again!')
         err=2;
      end
   end
end
if ~err
   if unitpick==16
      switch unit2
         % Frequency / Angular Velocity
      case 'GHz'
         x = x./1e9
      case 'MHz'
         x = x./1e6
      case 'kHz'
         x = x./1e3
      case 'Hz'
         x = x;
      case '1/min'
         x = x.*60;
      case '1/h'
         x = x.*3600;
      case 'rad/s'
         x = x.*(2*pi);
      case 'deg/s'
         x = x.*360;
      case 'rpm'
         x = x.*60;
      otherwise
         disp('Inconsistent unit - please check units and try again!')
         err=2;
      end
   end
end
if ~err
   if unitpick==17
      switch unit2
         % Angle
      case 'rad'
         x = x;
      case 'deg'
         x = x.*(180/pi);
      case 'rev'
         x = x./(2*pi);
      otherwise
         disp('Inconsistent unit - please check units and try again!')
         err=2;
      end
   end
end
if ~err
   if unitpick==18
      switch unit2
         % Stress Intensity Factor
      case 'MPa*m^1/2'
         x = x;
      case 'MPa*mm^1/2'
         x = x.*sqrt(1000);
      case 'MPa*in^1/2'
         x = x./6.894757./sqrt(0.0254);
      otherwise
         disp('Inconsistent unit - please check units and try again!')
         err=2;
      end
   end
end
if ~err
   if unitpick==19
      switch unit2
         % Fuel Consumption
      case 'l/100km'
         x = x;
      case 'miles/galUS'
         x = 235.214596754951./x;
      otherwise
         disp('Inconsistent unit - please check units and try again!')
         err=2;
      end
   end
end
if ~err
   if unitpick==20
      switch unit2
         % Viscosity
      case 'Ns/m^2'
         x = x;
      case 'poise'
         x = x.*10;
      case 'centipoise'
         x = x.*1000;
      case 'lbfs/sq-ft'
         x = x.*0.02089;
      otherwise
         disp('Inconsistent unit - please check units and try again!')
         err=2;
      end
   end
end
if ~err
   if unitpick==21
      switch unit2
         % Strain
      case 'm/m'
         x = x;
      case 'mm/m'
         x = x.*1e3;
      case 'ym/m'
         x = x.*1e6;
      case 'nm/m'
         x = x.*1e9;
      case '%'
         x = x.*1e2;
      case 'o/oo'
         x = x.*1e3;
      case 'in/in'
         x = x;
      case 'mill/in'
         x = x.*1e3';
      otherwise
         disp('Inconsistent unit - please check units and try again!')
         err=2;
      end
   end
end
if ~err
   if unitpick==22
      switch unit2
         % Metric Prefixes
      case 'Yotta'
         x = x.*1e-24;
      case 'Zetta'
         x = x.*1e-21;
      case 'Exa'
         x = x.*1e-18;
      case 'Peta'
         x = x.*1e-15;
      case 'Tera'
         x = x.*1e-12;
      case 'Giga'
         x = x.*1e-9;
      case 'Mega'
         x = x.*1e-6;
      case 'Myria'
         x = x.*1e-5;
      case 'kilo'
         x = x.*1e-3;
      case 'hecto'
         x = x.*1e-2;
      case '1'
         x = x;
      case 'deci'
         x = x.*1e1;
      case 'centi'
         x = x.*1e2;
      case 'milli'
         x = x.*1e3;
      case 'micro'
         x = x.*1e6;
      case 'nano'
         x = x.*1e9;
      case 'pico'
         x = x.*1e12;
      case 'femto'
         x = x.*1e15;
      case 'atto'
         x = x.*1e18;
      case 'zepto'
         x = x.*1e21;
      case 'yocto'
         x = x.*1e24;
      end
   end
end

if ~err
   answer=x;
end


