function bodyFunction = fuselage_coordinate_generator(type)

if nargin==0
    clc
    clear all
    bodyFunction = fuselage_coordinate_generator('cylinder');
    coord = bodyFunction(2.5,0.2);
    return
end

switch type
    case 'cylinder'
        bodyFunction = @cylinder;
end

    function coord = cylinder(length,diameter)
        R = diameter/2;
        cLen = length-diameter;
        nRad = 2;
        nCyl = 10;
        theta1 = 0:pi/(2*nRad):pi/2;
        theta1 = theta1(1:end-1);
        theta2 = pi/2:pi/(2*nRad):3/2*pi;
        theta2 = theta2(2:end-1);
        theta3 = -pi/2:pi/(2*nRad):0;
        theta3 = theta3(2:end);
        circ1 = [cos(theta1); sin(theta1)]*R;
        circ2 = [cos(theta2); sin(theta2)]*R;
        circ3 = [cos(theta3); sin(theta3)]*R;
        circ1(1,:) = circ1(1,:)+R+cLen;
        circ2(1,:) = circ2(1,:)+R;
        circ3(1,:) = circ3(1,:)+R+cLen;
        upLineX = 0:cLen/nCyl:cLen;
        upLineY =  R+upLineX*0;
        upLine = [upLineX+R; upLineY];
        loLine = [upLineX+R; -upLineY];
        coord = [circ1, fliplr(upLine), circ2, loLine, circ3]';
    end

end