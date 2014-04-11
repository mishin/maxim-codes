function generate_datcom_input_config()

config = load_configuration();
fc = get_datacom_flight_conditions();
fid = fopen('for005.dat','wt');

fprintf(fid,'CASEID 1 GBU-38 BETA=0.0,PHI=0\n');
fprintf(fid,'DIM M\n');
% Now these parameters are fixed.
fprintf(fid,'$FLTCON NALPHA=15.,\n');
fprintf(fid,' ALPHA=-8.,-4.,0.,4.,8.,12.,16.,20.,24.,28.,32.,\n'); 
fprintf(fid,' ALPHA(12)=36.,40.,44.,48.,\n');
fprintf(fid,' BETA=0.0,\n'); 
fprintf(fid,' NMACH=1.0,\n'); 
fprintf(fid,' MACH=0.8,\n'); 
fprintf(fid,' ALT=13700.0,$\n'); 
% fprintf(fid,' REN=300000.0,$\n')

fprintf(fid,'$REFQ\n');
fprintf(fid,' XCG=%.4f,\n',config.CG(1));
fprintf(fid,' LREF=%.4f,\n',config.body.length);
fprintf(fid,' SREF=%.8f,$\n',config.body.sref);

fprintf(fid,'$AXIBOD TNOSE=CONICAL,\n');
fprintf(fid,' LNOSE=%.4f,\n',config.body.length_nose);
fprintf(fid,' DNOSE=%.4f,\n',config.body.diameter_nose);
fprintf(fid,' BNOSE=%.4f,\n',config.body.bluntness_nose);
fprintf(fid,' LCENTR=%.4f,\n',config.body.length_center);
fprintf(fid,' DCENTR=%.4f,\n',config.body.diameter_center);
% fprintf(fid,' LAFT=%.4f,\n',config.body.length_aft);
% fprintf(fid,' DAFT=%.4f,\n',config.body.diameter_aft);
fprintf(fid,' DEXIT=%.4f,$\n',config.body.diameter_exit);

fprintf(fid, '$PROTUB NPROT=0.0,$\n');

fprintf(fid,'$FINSET1\n');
fprintf(fid,' NPANEL=%.1f,\n',config.wing.numberOfTails);
fprintf(fid,' ZUPPER=%.7f,\n',config.wing.zupper);
fprintf(fid,' ZLOWER=%.7f,\n',config.wing.zlower);
fprintf(fid,' SSPAN=%.1f,%.6f,\n',config.wing.span(1),config.wing.span(2));
fprintf(fid,' CHORD=%.6f,%.6f,\n',config.wing.secChords(1),config.wing.secChords(2));
fprintf(fid,' XLE=%.5f,\n',config.wing.xleading_edge_location);
fprintf(fid,' SWEEP=%.1f\n',config.wing.sweepAngle);
fprintf(fid,' STA=1.0,\n'); %If STA = 0, the sweep angle input is measured
PHIF = config.wing.PHIF;
n = length(PHIF);
fprintf(fid,' PHIF=');
for i=1:n-1
    fprintf(fid,'%.1f,',PHIF(i));
end
fprintf(fid,'%.1f,\n',PHIF(n));
fprintf(fid,'$END\n');

fprintf(fid,'$FINSET2\n');
fprintf(fid,' NPANEL=%.1f,\n',config.fin.numberOfTails);
fprintf(fid,' ZUPPER=%.7f,\n',config.fin.zupper);
fprintf(fid,' ZLOWER=%.7f,\n',config.fin.lower);
fprintf(fid,' LMAXU=%.7f,\n',config.fin.lmaxu);
fprintf(fid,' LMAXL=%.7f,\n',config.fin.lmaxl);
fprintf(fid,' LFLATU=%.6f,\n',config.fin.lflatu);
fprintf(fid,' LFLATL=%.6f,\n',config.fin.lflatl);
fprintf(fid,' SSPAN=%.3f,%.5f,\n',config.fin.span(1),config.fin.span(2));
fprintf(fid,' CHORD=%.6f,%.4f,\n',config.fin.secChords(1),config.fin.secChords(2));
fprintf(fid,' XLE=%.5f,\n',config.fin.xleading_edge_location);
fprintf(fid,' SWEEP=%.1f,\n',config.fin.sweepAngle);
fprintf(fid,' STA=1.0,\n');
PHIF = config.fin.PHIF;
n = length(PHIF);
fprintf(fid,' PHIF=');
for i=1:n-1
    fprintf(fid,'%.1f,',PHIF(i));
end
fprintf(fid,'%.1f,\n',PHIF(n));
fprintf(fid,'$END\n');

% fprintf(fid,'# ========== DEFLECTION ========== #\n');
delta = fc.deflection;
fprintf(fid,'$DEFLCT DELTA2=%.1f,%.1f,%.1f,%.1f,$\n',delta(1),delta(2),delta(3),delta(4));
% fprintf(fid,'PRINT GEOM BODY\n');
% fprintf(fid,'PRINT AERO BODY\n');
fprintf(fid,'DERIV DEG\n');
fprintf(fid,'PART\n');
fprintf(fid,'PLOT\n');
fprintf(fid,'DAMP\n');
fprintf(fid,'SOSE\n');
fprintf(fid,'SAVE\n');
fprintf(fid,'NEXT CASE\n');
fprintf(fid,'CASEID TRIM OF CASE NUMBER 1\n');
fprintf(fid,'$TRIM SET=2.0,$\n');
fprintf(fid,'PRINT AERO TRIM\n');
fprintf(fid,'NEXT CASE\n');

fclose(fid);

end
