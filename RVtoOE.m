function [a,e,i,BOmega,lomega,f,h,e_hat,n_hat,h_hat,t]  = RVtoOE(r_vec,v_vec,mu)
%this will take the state vectors and return the orbital elements and the
%cooresponding vectors 
% a is semi major axis
% e is eccentricity
% i is inclination
% BOmega is longitude of the ascending node
% lomega is argument of periapsis 
% f is true anomoly 
% mu is graviational parameter 

%equations for orbital elements
r = norm(r_vec);
r_hat = r_vec/r;

h_vec = cross(r_vec,v_vec);
h = norm(h_vec);
h_hat = h_vec/h;

e_vec = (1/mu)*cross(v_vec,h_vec)-r_hat;
e = norm(e_vec);
e_hat = e_vec/e;

n_hat = cross([0 0 1],h_vec)/norm(cross([0 0 1],h_vec));

epsilon = 0.5*dot(v_vec,v_vec)-(mu/r);

a = -mu/(2*epsilon);

i = acos(dot(h_hat,[0 0 1]));

BOmega = atan2(dot(n_hat,[0 1 0]),dot(n_hat,[1 0 0]));

if dot(e_hat,[0 0 1]) >= 0
    lomega = acos(dot(n_hat,e_hat));   
else
    lomega = 2*pi - acos(dot(n_hat,e_hat));    
end

if dot(r_vec,v_vec) >= 0   
    f = acos(dot(e_hat,r_hat));    
else
     f = 2*pi - acos(dot(e_hat,r_hat));
end

%finds the time
if e < 1
    E = 2*atan2(sqrt((1-e)/(1+e))*tan(f/2),1);
    M = E-e*sin(E);
    t = sqrt(a^3/mu)*M;
    F = NaN;
else
    F = 2*atanh(tan(f/2)*sqrt((e-1)/(e+1)));
    M = e*sinh(F)-F;
    t = sqrt(-a^3/mu)*M;
    E = NaN;
end

end

