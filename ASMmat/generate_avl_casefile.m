function generate_avl_casefile(files,flightConditions,setup)

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
fprintf(fid,'a a %.6f\n',setup.alpha);
fprintf(fid,'b b %.6f\n',setup.beta);
fprintf(fid,'d1 d1 %.6f\n',setup.fin1);
fprintf(fid,'d2 d2 %.6f\n',setup.fin2);
fprintf(fid,'d3 d3 %.6f\n',setup.fin3);
fprintf(fid,'d4 d4 %.6f\n',setup.fin4);
fprintf(fid,'x\n');
fprintf(fid,'st\n\n');
fprintf(fid,'st\n%s\n',stabilityOutput);
fprintf(fid,'\nquit\n');
fclose(fid);
end