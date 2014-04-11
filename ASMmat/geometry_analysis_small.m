function config = geometry_analysis_small(config)

    function wing = get_PHIF(wing)
        step = 360/wing.numberOfTails;
        i = 1;
        for angle = 0:step:360-step
            wing.PHIF(i) = angle + wing.xAngle;
            i = i+1;
        end
    end

config.wing = get_PHIF(config.wing);
config.fin = get_PHIF(config.fin);

end