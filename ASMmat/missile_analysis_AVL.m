function missile_analysis_AVL(missile,flightConditions)
clc
% 1. generate input file
% 2. run avl
% 3. collect output
% 4. process output

if nargin==0
    fc.Mach = 0.3;
    fc.velocity = 100;
    fc.density = 1.05;
    missile = geometry_analysis();
    missile_analysis_AVL(missile,fc);
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

%dosCommand=sprintf('%s\\avl.exe %s < %s ',avlPath,files.input,files.case);
dosCommand=sprintf('%s %s < %s ',avlPath,files.input,files.case);
[~,~] = dos(dosCommand);
results = collect_stability_results(files.output);

delete(files.input);
delete(files.case);
%delete(files.output);
end