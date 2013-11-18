function result = collect_stability_results(filePath)

if nargin==0
    collect_stability_results('xtail.stab');
    return
end

out = fileread(filePath);

iAlpha    =strfind(out,'Alpha')+7;
iElevator =strfind(out,'elevator')+17;
iNP       =strfind(out,'Neutral')+20;
ie        =strfind(out,'e =')+3;
iCL       =strfind(out,'CLtot')+7;
iCD       =strfind(out,'CDtot')+7;
iCl       =strfind(out,'Cltot')+7;
iCn       =strfind(out,'Cntot')+7;
iCm       =strfind(out,'Cmtot')+7;
iCDind    =strfind(out,'CDind')+7;
iSpiral   =strfind(out,'Clb Cnr / Clr Cnb')+20;

iCLa      =strfind(out,'CLa')+5;
iCYa      =strfind(out,'CYa')+5;
iCla      =strfind(out,'Cla')+5;
iCma      =strfind(out,'Cma')+5;
iCna      =strfind(out,'Cna')+5;
iCLb      =strfind(out,'CLb')+5;
iCYb      =strfind(out,'CYb')+5;
iClb      =strfind(out,'Clb')+5;
iCmb      =strfind(out,'Cmb')+5;
iCnb      =strfind(out,'Cnb')+5;

iCLp      =strfind(out,'CLp')+5;
iCYp      =strfind(out,'CYp')+5;
iClp      =strfind(out,'Clp')+5;
iCmp      =strfind(out,'Cmp')+5;
iCnp      =strfind(out,'Cnp')+5;
iCLq      =strfind(out,'CLq')+5;
iCYq      =strfind(out,'CYq')+5;
iCmq      =strfind(out,'Cmq')+5;
iCnq      =strfind(out,'Cnq')+5;
iCLr      =strfind(out,'CLr')+5;
iCYr      =strfind(out,'CYr')+5;
iClr      =strfind(out,'Clr')+5;
iCmr      =strfind(out,'Cmr')+5;
iCnr      =strfind(out,'Cnr')+5;

iCLd1     =strfind(out,'CLd1')+6;
iCYd1     =strfind(out,'CYd1')+6;
iCld1     =strfind(out,'Cld1')+6;
iCmd1     =strfind(out,'Cmd1')+6;
iCnd1     =strfind(out,'Cnd1')+6;
iCDffd1   =strfind(out,'CDffd1')+8;
ied1      =strfind(out,'ed1')+5;

iCLd2     =strfind(out,'CLd2')+6;
iCYd2     =strfind(out,'CYd2')+6;
iCld2     =strfind(out,'Cld2')+6;
iCmd2     =strfind(out,'Cmd2')+6;
iCnd2     =strfind(out,'Cnd2')+6;
iCDffd2   =strfind(out,'CDffd2')+8;
ied2      =strfind(out,'ed2')+5;

iCLd3     =strfind(out,'CLd3')+6;
iCYd3     =strfind(out,'CYd3')+6;
iCld3     =strfind(out,'Cld3')+6;
iCmd3     =strfind(out,'Cmd3')+6;
iCnd3     =strfind(out,'Cnd3')+6;
iCDffd3   =strfind(out,'CDffd3')+8;
ied3      =strfind(out,'ed3')+5;

iCLd4     =strfind(out,'CLd4')+6;
iCYd4     =strfind(out,'CYd4')+6;
iCld4     =strfind(out,'Cld4')+6;
iCmd4     =strfind(out,'Cmd4')+6;
iCnd4     =strfind(out,'Cnd4')+6;
iCDffd4   =strfind(out,'CDffd4')+8;
ied4      =strfind(out,'ed4')+5;

result.alpha = str2double(out(iAlpha:iAlpha+10))*pi/180;

end