function [deltav,TOF,at1,at2] = delv(r1,r2,r3,del,mu,trigger)
    
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
    
elseif trigger == 'phs'
    
    n = sqrt(mu/r1^3);
    TOF = (2*pi+del)/n;
    at1 = (mu*(TOF/(2*pi))^2)^(1/3);
    deltav = 2*abs(-sqrt(mu*(2/r1-1/at1))+sqrt(mu/r1));
    
    at2 = 0;

else 
    
    fprintf('what\n');
      
end

end

