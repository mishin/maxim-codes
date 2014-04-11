function AVL_batch_run()
clc
clear all

Mach = [0.1, 0.25, 0.4, 0.55, 0.7];
alpha = [-15 -12 -9 -6 -3 0 3 6 9 12 15];
elevator = [-30 -20 -10 0 10 20 30];
altitude = 0;

savefile = 'AVL_analysis_20131129.txt';

n1 = length(Mach);
n2 = length(alpha);
n3 = length(elevator);
ntot = n1*n2*n3;
designMatrix = zeros(ntot,3);

%% design matrix
curNumber = 1;
for i=1:n1
    for j=1:n2
        for k=1:n3
            designMatrix(curNumber,1) = Mach(i);
            designMatrix(curNumber,2) = alpha(j);
            designMatrix(curNumber,3) = elevator(k);
            curNumber = curNumber +1;
        end
    end
end

%% aerodynamic analysis
CL = zeros(ntot,1);
CD = zeros(ntot,1);
CM = zeros(ntot,1);
control = control_deflections;
gb = load_configuration(1);
gb.CG = [0.8148, 0.0, 0.0];
fid = fopen(savefile,'wt');
fprintf(fid,'Mach\tVelocity\talpha\televator\tCL\tCD\tCM\n');
fclose(fid);
for i=1:ntot
    fc = get_flight_conditions(altitude,designMatrix(i,1));
    setup.alpha = designMatrix(i,2);
    setup.beta = 0;
    setup = control.set_de(setup,designMatrix(i,3));
    aero = missile_analysis_AVL(gb,fc,setup);
    CL(i) = aero.CL;
    CD(i) = aero.CD;
    CM(i) = aero.Cm;
    fid = fopen(savefile,'a');
    fprintf(fid,'%.2f\t%.2f\t%.1f\t%.1f\t%.6e\t%.6e\t%.6e\n',fc.Mach,fc.velocity,designMatrix(i,2),designMatrix(i,3),CL(i),CD(i),CM(i));
    fclose(fid);
end

end