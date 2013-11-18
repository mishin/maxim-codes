function generate_avl_casefile(files,flightConditions)

filePath = files.case;
stabilityOutput = files.output;

fid = fopen(filePath,'wt');
fprintf(fid,'OPER\n');
%fprintf(fid,'g\n\n');
fprintf(fid,'m\n');
fprintf(fid,'mn %.4f\n',flightConditions.Mach);
fprintf(fid,'v %.4f\n',flightConditions.velocity);
fprintf(fid,'g 9.80665\n'); % gravity acceleration is defined here, this statement should be proved!!!
fprintf(fid,'\n');
fprintf(fid,'x\n');
fprintf(fid,'st\n\n');
fprintf(fid,'st\n%s\n',stabilityOutput);
fprintf(fid,'\nquit\n');
fclose(fid);
end