function control = control_deflections()

control.defaults = @get_defalults;
control.set_de = @set_elevator_deflection;
control.set_da = @set_aileron_deflection;
control.set_dr = @set_rudder_deflection;

    function control = get_defalults()
        control.alpha = 0;
        control.beta = 0;
        control.fin1 = 0;
        control.fin2 = 0;
        control.fin3 = 0;
        control.fin4 = 0;
    end
    function control = set_elevator_deflection(control,deflection)
        control.fin1 = deflection;
        control.fin2 = -deflection;
        control.fin3 = -deflection;
        control.fin4 = deflection;
    end

    function control = set_aileron_deflection(control,deflection)
        control.fin1 = -deflection;
        control.fin2 = -deflection;
        control.fin3 = -deflection;
        control.fin4 = -deflection;
    end

    function control = set_rudder_deflection(control,deflection)
        control.fin1 = deflection;
        control.fin2 = deflection;
        control.fin3 = -deflection;
        control.fin4 = -deflection;
    end
end