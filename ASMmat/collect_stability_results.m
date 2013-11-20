function result = collect_stability_results(filePath)

if nargin==0
    collect_stability_results('xtail.stab');
    return
end

out = fileread(filePath);
radDeg = rad2deg(1);

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
iCDff     =strfind(out,'CDff')+7;
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

result.alpha        = str2double(out(iAlpha:iAlpha+10));%*pi/180;
result.elevator     = str2double(out(iElevator:iElevator+10))*pi/180;
result.xNP          = str2double(out(iNP:iNP+11));
result.CL           = str2double(out(iCL:iCL+10));
result.CD           = str2double(out(iCD:iCD+10));
result.CDind        = str2double(out(iCDind:iCDind+10));
result.CDff         = str2double(out(iCDff:iCDff+10));
result.k            = result.CDff/result.CL^2;

result.e         = str2double(out(ie:ie+10));            
result.Cl        = str2double(out(iCl:iCl+10));
result.Cn        = str2double(out(iCn:iCn+10));
result.Cm        = str2double(out(iCm:iCm+10));
result.spiral    = str2double(out(iSpiral:iSpiral+11));
result.derivs.CLa = str2double(out(iCLa:iCLa+11)) /radDeg;
result.a          = result.derivs.CLa;
result.derivs.CYa = str2double(out(iCYa:iCYa+11)) /radDeg;
result.derivs.Cla = str2double(out(iCla:iCla+11)) /radDeg;
result.derivs.Cma = str2double(out(iCma:iCma+11)) /radDeg;
result.derivs.Cna = str2double(out(iCna:iCna+11)) /radDeg;
result.derivs.CLb = str2double(out(iCLb:iCLb+11)) /radDeg;
result.derivs.CYb = str2double(out(iCYb:iCYb+11)) /radDeg;
result.derivs.Clb = str2double(out(iClb:iClb+11)) /radDeg;
result.derivs.Cmb = str2double(out(iCmb:iCmb+11)) /radDeg;
result.derivs.Cnb = str2double(out(iCnb:iCnb+11)) /radDeg;
result.CL0 = result.CL-result.derivs.CLa*result.alpha;

result.derivs.CLp=str2double(out(iCLp:iCLp+11));
result.derivs.CYp=str2double(out(iCYp:iCYp+11));
result.derivs.Clp=str2double(out(iClp:iClp+11));
result.derivs.Cmp=str2double(out(iCmp:iCmp+11));
result.derivs.CLp=str2double(out(iCnp:iCnp+11));
result.derivs.CLq=str2double(out(iCLq:iCLq+11));
result.derivs.CYq=str2double(out(iCYq:iCYq+11));
result.derivs.Cmq=str2double(out(iCmq:iCmq+11));
result.derivs.Cnq=str2double(out(iCnq:iCnq+11));
result.derivs.CLr=str2double(out(iCLr:iCLr+11));
result.derivs.CYr=str2double(out(iCYr:iCYr+11));
result.derivs.Clr=str2double(out(iClr:iClr+11));
result.derivs.Cmr=str2double(out(iCmr:iCmr+11));
result.derivs.Cnr=str2double(out(iCnr:iCnr+11));

result.control.CLd1=str2double(out(iCLd1:iCLd1+11));
result.control.CYd1=str2double(out(iCYd1:iCYd1+11));
result.control.Cld1=str2double(out(iCld1:iCld1+11));
result.control.Cmd1=str2double(out(iCmd1:iCmd1+11));
result.control.Cnd1=str2double(out(iCnd1:iCnd1+11));
result.control.CLd1=str2double(out(iCDffd1:iCDffd1+11));
result.control.ed1=str2double(out(ied1:ied1+11));

result.control.CLd2=str2double(out(iCLd2:iCLd2+11));
result.control.CYd2=str2double(out(iCYd2:iCYd2+11));
result.control.Cld2=str2double(out(iCld2:iCld2+11));
result.control.Cmd2=str2double(out(iCmd2:iCmd2+11));
result.control.Cnd2=str2double(out(iCnd2:iCnd2+11));
result.control.CLd2=str2double(out(iCDffd2:iCDffd2+11));
result.control.ed2=str2double(out(ied2:ied2+11));

result.control.CLd3=str2double(out(iCLd3:iCLd3+11));
result.control.CYd3=str2double(out(iCYd3:iCYd3+11));
result.control.Cld3=str2double(out(iCld3:iCld3+11));
result.control.Cmd3=str2double(out(iCmd3:iCmd3+11));
result.control.Cnd3=str2double(out(iCnd3:iCnd3+11));
result.control.CLd3=str2double(out(iCDffd3:iCDffd3+11));
result.control.ed3=str2double(out(ied3:ied3+11));

result.control.CLd4=str2double(out(iCLd4:iCLd4+11));
result.control.CYd4=str2double(out(iCYd4:iCYd4+11));
result.control.Cld4=str2double(out(iCld4:iCld4+11));
result.control.Cmd4=str2double(out(iCmd4:iCmd4+11));
result.control.Cnd4=str2double(out(iCnd4:iCnd4+11));
result.control.CLd4=str2double(out(iCDffd4:iCDffd4+11));
result.control.ed4=str2double(out(ied4:ied4+11));

result.control.CLde = result.control.CLd1 - result.control.CLd2 + result.control.CLd3 - result.control.CLd4;


end