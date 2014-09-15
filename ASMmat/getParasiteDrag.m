function result=getParasiteDrag(aircraft,altitude,mach)
if exist('aircraft','var')~=1
    result=getParasiteDragTest();
    return
end

%option=2;
%pathDat=load('LSApathData.mat');pathDat=pathDat.LSApathData;
%aero_Cd0_path=pathDat.aero_CD0;
%inputFile ='PDragIn.inp';
%outputFile='PDragOut.out';
inputFile = 'drag_p_ktx2.inp';
outputFile = 'drag_p_ktx2.out';
%% Test function
% Read inputs for friction drag calculation
M         =mach;
h         =altitude;
W         =aircraft.weight; % kg  
[rho,a]   =ISAtmosphere(h); 
Sref      =aircraft.wing.area;
CL=W/(0.5*rho*(M*a)^2*Sref);

%% Inputs for Cdo Calculations (Note: Unit in Ft)  
a1=3   ;                    % No of Components (Wing, Fuselage, HT) 
a2=units(altitude,'m','ft');% Flight Altitude in ft 
a3=units(aircraft.wing.area,'m^2','sq-ft'); % Wing area (in ft^2) 
a4=0.095;                   % Base area  ??????????????????????????????
a5=0; % nacelle exit are
a6=1.8; % Experiment factor,(Sears-Haack body? ?? ??? ?? ??? ?)
a7=10;% Leakage drag
%b1=aircraft.wing.LEsweep*180/pi; %Leading edge sweepback angle
b1 = aircraft.wing.sweepLE;
b2=aircraft.wing.sweepC4;
b3 = aircraft.wing.sweepC2;
%b2=aircraft.wing.C4sweep*180/pi;
%b3=aircraft.wing.C2sweep*180/pi;
b4=units(aircraft.wing.span,'m','ft');% Note 
b5=pi*units(aircraft.body.diameter,'m','ft')^2/4;% body max. area
b6=aircraft.body.diameter/aircraft.body.length; % Thickness ratio
%b7=aircraft.wing.maxCamber; %max. camber
b7 = 0.0;
b8=CL; %section lift coefficient
b9=1.02; % Constant 
% Wing input
w1=units(aircraft.wing.wettedArea,'m^2','sq-ft');
w2=units(aircraft.wing.span,'m','ft');%ref(i)	42.2	Reference length (ft)
w3=1;%stype(i)	2.	form factor type
w4=aircraft.wing.tc; %tc(i)	8.45	component length/sqrt(width*height)
w5=aircraft.wing.sweepC4;%*180/pi;%amt(i)	0.	?? ??????? ??? (sweep angle at max thickness ratio 
w6=aircraft.wing.tc;%xt(i)	0.	?? ???? ?? (max Thickness ratio)
w7=1; %Q(i)	1.	

% Fuselage input 
f1=units(aircraft.body.wettedArea,'m^2','sq-ft');
f2=units(aircraft.body.length,'m','ft');
f3=2 ;
f4=aircraft.body.length/aircraft.body.diameter;
f5=0;
f6=0;
f7=1;

%HT and VTail inputs
HT1=units(2*aircraft.fin.wettedArea,'m^2','sq-ft');
HT2=units(2*aircraft.fin.equivSpan,'m','ft');
HT3=1;
HT4=aircraft.fin.thickness;
HT5=aircraft.fin.sweepC4;
HT6=aircraft.fin.thickness;
HT7=1.08;

%% Run Aero execuation file from Fortran 
%if option==1
    fID=fopen(inputFile,'w');
    fprintf(fID,'%.1f %.1f %.1f %.5f %.1f %.1f %.1f\n',a1,a2,a3,a4,a5,a6,a7);
    fprintf(fID,'%.1f %.1f %.1f %.1f %.2f %.2f %.5f %.1f %.2f\n',b1,b2,b3,b4,b5,b6,b7,b8,b9);
    fprintf(fID,'%s %.2f %.2f %.2f %.2f %.2f %.2f %.2f\n','Wing',w1,w2,w3,w4,w5,w6,w7);
    fprintf(fID,'%s %.2f %.2f %.2f %.2f %.2f %.2f %.2f\n','Fuselage',f1,f2,f3,f4,f5,f6,f7);
    fprintf(fID,'%s %.3f %.3f %.3f %.3f %.3f %.3f %.3f\n','H-Tail',HT1,HT2,HT3,HT4,HT5,HT6,HT7);
    fclose(fID);
    lastPath=pwd;
    %cd(aero_Cd0_path)
    %dos('aero_CD0.exe');
    dos('Min_Drag.exe');
    %cd(lastPath);
    % Read outputs from AERO 
    fidVariable=fopen(outputFile);
    data=textscan(fidVariable,'%f %f %f %f %f %f %f %f','HeaderLines',3);
    fclose(fidVariable);
    Mref=data{1}; 
    CDoRef=data{8};
    %CDo2=spline(Mref(1:end-1),CDoRef(1:end-1),M);
    CDo2 = interp1(Mref(1:end-1), CDoRef(1:end-1),M,'pchip');
    delete(inputFile, outputFile);
    %figure(3)
    %plot(Mref, CDoRef)
%else
%    CDo2=getCd0(aircraft,mach,altitude);
%end

result   =CDo2;

%% Test Function
    function result=getParasiteDragTest()
        clc
        close all
        aircraft=load_configuration_large();
        mArray=.2:.01:1.0;
        i=0;
        for Mach=mArray
            i=i+1;
            Cd0(i)=getParasiteDrag(aircraft,1500,Mach);
        end
        plot(mArray,Cd0)
        result=Cd0(end);
    end
end
