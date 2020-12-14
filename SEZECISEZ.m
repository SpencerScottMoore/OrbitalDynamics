function [conout] = SEZECISEZ(lambda,phi,convert,trigger)
    %defines the rotation matrix and inverse. Origininally just used matlab
    %built in inv() function but typed it out during trouble shooting phase
    %so i just left it
    ROT = [sin(phi) 0 -cos(phi); 0 1 0; cos(phi) 0 sin(phi)]*[cos(lambda)...
        sin(lambda) 0; -sin(lambda) cos(lambda) 0; 0 0 1];
    ROTin = [cos(lambda) -sin(lambda) 0; sin(lambda) cos(lambda) 0; 0 0 1]...
        *[sin(phi) 0 cos(phi); 0 1 0; -cos(phi) 0 sin(phi)];

    %figures out whether user wanted to convert from SEZ to ECI or vice
    %versa by the arguments specifiying desired output
    if trigger == 'SEZ'
        conout = ROT*convert;
    end
    
    if trigger == 'ECI'
        conout = ROTin*convert;
    end
end

