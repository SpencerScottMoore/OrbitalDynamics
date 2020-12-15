function [f,E,F] = NewtonMethod(M,tol,e)
%This function uses newtons method to determine the time of flight from the
%mean anamoly, a accuracy tolerance (tol = 0.01 works well), and orbits 
%eccentricity. the code will output 
%f, true anomaly
%E, eccentric anomoly, if applicable
%F, hyperbolic anomoly, if applicable


if e < 1
    %initial guess
    E0 = M;
    E(1)= E0;
    error = 1;
    p=1;

    %goes until value is close
    while error > tol

      E(p+1) = E(p)-(E(p)-e*sin(E(p))-M)/(1-e*cos(E(p)));   
      error = abs(E(p)-E(p+1));
      p=p+1;

    end

    %outputs f
    f = 2*atan2(sqrt((1+e)/(1-e))*tan(E(end)/2),1);
    F = NaN;
else
    %initial guess
    F0 = M;
    F(1)= F0;
    error = 1;
    p=1;

    %goes until value is close
    while error > tol

      F(p+1) = F(p)-(-F(p)+e*sinh(F(p))-M)/(e*cosh(F(p))-1);   
      error = abs(F(p)-F(p+1));
      p=p+1;

    end

    %outputs f
    f = 2*atan(sqrt((e+1)/(e-1))*tanh(F(end)/2));
    E=NaN;
end

end

