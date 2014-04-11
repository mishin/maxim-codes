function massOutput = weight_analysis_large(gb,massInput)

if nargin==0
   gb = load_configuration(1);
   mIn = get_default_mass_input();
   mIn.Vmax = gb.Vmax;
   massOutput = weight_analysis_large(gb,mIn);
   return
elseif nargin==1
   mIn = get_default_mass_input();
   mIn.Vmax = gb.Vmax;
   massOutput = weight_analysis_large(gb,mIn);
   return
end

    function MassIn = get_default_mass_input()
        MassIn.payloadFusRatio = 0.5;           % x-payload / fuselage length
        MassIn.fuselageCGratio = 0.45;          % fuselage CGx / fuselage length
        MassIn.Xcg_seeker = 0.04425;            % x-location of seeker (m)   
        MassIn.Xcg_imu = 0.3135;                % x-location of IMU (m)
        MassIn.Xcg_gps = 0.3135;                % x-location of GPS (m)
        MassIn.Xcg_seekele = 0.5385;            % x-location of seeker electronics (m)
        MassIn.Xcg_misscom = 1.35624;           % x-location of mission computer (m)
        MassIn.Xcg_batt = 1.47882;              % x-location of thermal battery (m)
        MassIn.Xcg_cas = 1.6899;                % x-location of control actuator system (m)
    end

massOutput = GBMass(gb,massInput);

end