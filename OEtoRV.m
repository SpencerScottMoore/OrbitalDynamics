function  [r_vec,v_vec]= OEtoRV(a,e,i,BOmega,lomega,f,mu)
%this function will take orbital elements and return the state vectors of
%the orbit
% a is semi major axis
% e is eccentricity
% i is inclination
% BOmega is longitude of the ascending node
% lomega is argument of periapsis 
% f is true anomoly 
% mu is graviational parameter 


%find h
h = sqrt(mu*a*(1-e^2));

%finds r and v in PQW frame
r_PQW_vec = (a*(1-e^2))/(1+e*cos(f))*[cos(f) ; sin(f); 0];
v_PQW_vec = mu/h*[-sin(f); e+cos(f); 0];

%Rotation matrices
PQWtoECIBO = [cos(-BOmega), sin(-BOmega), 0; -sin(-BOmega), cos(-BOmega), 0; 0, 0, 1];
PQWtoECIi = [1, 0, 0; 0, cos(-i), sin(-i); 0, -sin(-i), cos(-i)];
PQWtoECIlo = [cos(-lomega), sin(-lomega), 0; -sin(-lomega), cos(-lomega), 0; 0, 0, 1];

%convert to ECI frame rotation matrix
PQWtoECI = PQWtoECIBO*PQWtoECIi*PQWtoECIlo;

%Converts to final r and v
r_vec = PQWtoECI*r_PQW_vec;
v_vec = PQWtoECI*v_PQW_vec;

end

