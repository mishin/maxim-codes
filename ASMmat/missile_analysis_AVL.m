function results = missile_analysis_AVL(missile,flightConditions,setup)

if nargin==0
    clc
    alt = 5000;
    vel = 100;
    fc = get_flight_conditions(alt,vel);
    missile = load_configuration_large();
    control = control_deflections();
    setup.alpha = 0;
    setup.beta = 0;
    setup = control.set_de(setup,0);
    results = missile_analysis_AVL(missile,fc,setup);
    return
end

avlPath = 'avl_win32.exe';
name = 'xtail';
files.body = strcat(name,'.fus');
files.input = strcat(name,'.avl');
files.case = strcat(name,'.case');
files.output = strcat(name,'.stab');

missile.CD0 = getParasiteDrag(missile, flightConditions.altitude, flightConditions.Mach);
generate_avl_input_config(missile,files);
generate_avl_casefile(files,flightConditions,setup);

dosCommand=sprintf('%s %s < %s ',avlPath,files.input,files.case);
[~,~] = dos(dosCommand);
results = collect_stability_results(files.output,missile);
%fprintf('.')

delete(files.input);
delete(files.body);
delete(files.case);
delete(files.output);

end