function generate_avl_input_config(missile,filePath)

if nargin==0
    missile = geometry_analysis();
    path = 'xtail.avl';
    generate_avl_input_config(missile,path);
    return
end

body = missile.body;
wing = missile.wing;
fin = missile.fin;
panels = missile.VLMpanels;

fid = fopen(filePath.body,'wt');
for i=1:length(body.coord)
    fprintf(fid,'%.6f  %.6f\n',body.coord(i,1),body.coord(i,2));
end
fclose(fid);

fid = fopen(filePath.input, 'wt');
fprintf(fid,'Xtail analysis generated by AADL\n');
fprintf(fid,'0\n0 0 0\n');
fprintf(fid,'%.4f %.4f %.4f\n',wing.area, wing.MAC, wing.span);
fprintf(fid,'%.4f %.4f %.4f\n',missile.CG(1),missile.CG(2),missile.CG(3));
fprintf(fid,'%.6f\n',missile.CD0);

fprintf(fid,'# ========== BODY ==========\n');
fprintf(fid,'BODY\n');
fprintf(fid,'Fuselage\n');
fprintf(fid,'28 2\n');
fprintf(fid,'BFIL\n');
fprintf(fid,'%s\n',filePath.body);

fprintf(fid,'# ========== WING ==========\n');
fprintf(fid,'SURFACE\nwing\n');
fprintf(fid,'%d %d %d %d\n',panels.wing(1), 0, panels.wing(2), 0);
fprintf(fid,'YDUPLICATE\n0\n');
fprintf(fid,'ANGLE\n%.4f\n',0); % wing incidence
fprintf(fid,'TRANSLATE\n');
fprintf(fid,'%.4f %.4f %.4f\n',wing.location(1),0,wing.location(2));
fprintf(fid,'# ---------- wing section 1 ----------\n');
fprintf(fid,'SECTION\n');
fprintf(fid,'%.4f %.4f %.4f %.4f %.4f\n',wing.rootSection(1),wing.rootSection(2),wing.rootSection(3),wing.secChords(1),wing.secIncidence(1));
if wing.airfoil~=0
    fprintf(fid,'AFILE\n%s\n',wing.airfoil);
end
fprintf(fid,'# ---------- wing section 2 ----------\n');
fprintf(fid,'SECTION\n');
fprintf(fid,'%.4f %.4f %.4f %.4f %.4f\n',wing.tipSection(1),wing.tipSection(2),wing.tipSection(3),wing.secChords(2),wing.secIncidence(2));
if wing.airfoil~=0
    fprintf(fid,'AFILE\n%s\n',wing.airfoil);
end
write_tail(fid,fin,1);
write_tail(fid,fin,2);
write_tail(fid,fin,3);
write_tail(fid,fin,4);

    function write_tail(fileID,fin,n)
        fprintf(fileID,'# ========== TAIL %d ==========\n',n);
        fprintf(fileID,'SURFACE\n');
        fprintf(fileID,'Tail%d\n',n);
        fprintf(fileID,'%d %d %d %d\n',panels.tail(1),0,panels.tail(2),0);
        fprintf(fileID,'TRANSLATE\n');
        fprintf(fileID,'%.4f %.4f %.4f\n',fin.locationX,0,0);
        fprintf(fileID,'# ---------- tail section 1 ----------\n');
        fprintf(fileID,'SECTION\n');
        fprintf(fileID,'%.4f %.4f %.4f %.4f %.4f\n',fin.rootSection(n,1),fin.rootSection(n,2),fin.rootSection(n,3),fin.secChords(1),0);
        if fin.airfoil~=0
            fprintf(fileID,'AFILE\n%s\n',fin.airfoil);
        end
        fprintf(fileID,'# ---------- tail section 2 ----------\n');
        fprintf(fileID,'SECTION\n');
        fprintf(fileID,'%.4f %.4f %.4f %.4f %.4f\n',fin.middleSection(n,1),fin.middleSection(n,2),fin.middleSection(n,3),fin.middleChord,0);
        if fin.airfoil~=0
            fprintf(fileID,'AFILE\n%s\n',fin.airfoil);
        end
        fprintf(fileID,'CONTROL\n');
        fprintf(fileID,'elevator%d ',n);

        fprintf(fileID,'1 %.4f 0 0 0 +1\n',fin.elevatorChordRatio);
        fprintf(fileID,'# ---------- tail section 3 ----------\n');
        fprintf(fileID,'SECTION\n');
        fprintf(fileID,'%.4f %.4f %.4f %.4f %.4f\n',fin.tipSection(n,1),fin.tipSection(n,2),fin.tipSection(n,3),fin.secChords(2),0);
        if fin.airfoil~=0
            fprintf(fileID,'AFILE\n%s\n',fin.airfoil);
        end
        fprintf(fileID,'CONTROL\n');
        fprintf(fileID,'elevator%d ',n);
        fprintf(fileID,'1 %.4f 0 0 0 +1\n',fin.elevatorChordRatio);
    end
fclose(fid);
end