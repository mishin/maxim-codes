function results = missile_analysis_AVL(missile,flightConditions)
clc
% 1. generate input file
% 2. run avl
% 3. collect output
% 4. process output

if nargin==0
    alt = 5000;
    vel = 100;
    fc = get_flight_conditions(alt,vel);
    missile = geometry_analysis();
    results = missile_analysis_AVL(missile,fc);
    return
end

avlPath = 'avl_win32.exe';
name = 'xtail';
files.input = strcat(name,'.avl');
files.case = strcat(name,'.case');
files.output = strcat(name,'.stab');
files.case
files.output
files.input
generate_avl_input_config(missile,files.input);
generate_avl_casefile(files,flightConditions);

dosCommand=sprintf('%s %s < %s ',avlPath,files.input,files.case);
[~,~] = dos(dosCommand);
results = collect_stability_results(files.output);

%delete(files.input);
%delete(files.case);
%delete(files.output);

end