function generate_avl_input_config(missile,filePath)

if nargin==0
    missile = geometry_analysis();
    path = 'xtail.avl';
    generate_avl_input_config(missile,path);
    return
end

body = missile.body;
wing = missile.wing;
tail = missile.tail;
panels = missile.VLMpanels;

fid = fopen(filePath, 'wt');
fprintf(fid,'Xtail analysis generated by AADL\n');
fprintf(fid,'0\n0 0 0\n');
fprintf(fid,'%.4f %.4f %.4f\n',wing.area, wing.MAC, wing.span);
fprintf(fid,'%.4f %.4f %.4f\n',missile.CG(1),missile.CG(2),missile.CG(3));
fprintf(fid,'%.6f\n',missile.CD0);

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
fprintf(fid,'# ---------- wing section 2 ----------\n');
fprintf(fid,'SECTION\n');
fprintf(fid,'%.4f %.4f %.4f %.4f %.4f\n',wing.tipSection(1),wing.tipSection(2),wing.tipSection(3),wing.secChords(2),wing.secIncidence(2));

write_tail(fid,tail,1);
write_tail(fid,tail,2);
write_tail(fid,tail,3);
write_tail(fid,tail,4);

    function write_tail(fileID,tail,n)
        fprintf(fileID,'# ========== TAIL %d ==========\n',n);
        fprintf(fileID,'SURFACE\n');
        fprintf(fileID,'Tail%d\n',n);
        fprintf(fileID,'%d %d %d %d\n',panels.tail(1),0,panels.tail(2),0);
        fprintf(fileID,'TRANSLATE\n');
        fprintf(fileID,'%.4f %.4f %.4f\n',tail.locationX,0,0);
        fprintf(fileID,'# ---------- tail section 1 ----------\n');
        fprintf(fileID,'SECTION\n');
        fprintf(fileID,'%.4f %.4f %.4f %.4f %.4f\n',tail.rootSection(n,1),tail.rootSection(n,2),tail.rootSection(n,3),tail.secChords(1),0);
        fprintf(fileID,'CONTROL\n');
        fprintf(fileID,'elevator%d ',n);
        fprintf(fileID,'1 %.4f 0 0 0 +1\n',tail.elevatorChordRatio);
        fprintf(fileID,'# ---------- tail section 2 ----------\n');
        fprintf(fileID,'SECTION\n');
        fprintf(fileID,'%.4f %.4f %.4f %.4f %.4f\n',tail.tipSection(n,1),tail.tipSection(n,2),tail.tipSection(n,3),tail.secChords(2),0);
        fprintf(fileID,'CONTROL\n');
        fprintf(fileID,'elevator%d ',n);
        fprintf(fileID,'1 %.4f 0 0 0 +1\n',tail.elevatorChordRatio);
    end
fclose(fid);
end