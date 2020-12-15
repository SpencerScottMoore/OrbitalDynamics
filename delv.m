function [deltav,TOF,at1,at2] = delv(r1,r2,r3,del,mu,trigger)
%This function caclulates the delta v cost of 3 different orbital manuevers
%It will calculate the cost for a homann transfer, bi-elliptic transfer, and
%a phasing orbit. If a varible is not used 
%r1 is the radius of the initial orbit of the craft for all cases
%
%r2 is the final radius of the orbit achieved by the hohmann transfer or 
%the max distance of the bielliptic transfer orbits
%
%r3 is the final orbit of a craft that used the bielliptic manuever
%
%del is used only for phasing orbit and it is the angle between the target,
%the orbiting planet, and the spacecraft. Positive if the target is ahead,
%negative if behind. Measured in radians
%
%mu is the gravitational parameter. For earth it is mu = 398600
%
%the trigger gives the function which manuever is being used. 'hoh' hohmann
%transfer, 'bie' for bielliptic transfer, 'pha' for phasing orbit


if trigger == 'hoh'
    at1 = (r1+r2)/2;
    
    deltav = abs(sqrt(mu*(2/r1-1/at1))-sqrt(mu/r1))...
        + abs(-sqrt(mu*(2/r2-1/at1))+sqrt(mu/r2));
    
    TOF = pi*sqrt(at1^3/mu);
    
    at2 = 0;
      
elseif trigger == 'bie'
        
     at1 = (r1+r3)/2;
     at2 = (r2+r3)/2;
        
    deltav = abs(sqrt(mu*(2/r1-1/at1))-sqrt(mu/r1))...
        + abs(sqrt(mu*(2/r3-1/at2))- sqrt(mu*(2/r3-1/at1)))...
        + abs(-sqrt(mu*(2/r2-1/at2))+sqrt(mu/r2));
    
    TOF = pi*(sqrt(at1^3/mu)+sqrt(at2^3/mu));
    
elseif trigger == 'pha'
    
    n = sqrt(mu/r1^3);
    TOF = (2*pi+del)/n;
    at1 = (mu*(TOF/(2*pi))^2)^(1/3);
    deltav = 2*abs(-sqrt(mu*(2/r1-1/at1))+sqrt(mu/r1));
    
    at2 = 0;

else 
    
    fprintf('trigger not recognized\n');
      
end

end

