function r = NozzleContour(dt,de,theta_n,theta_e,points)

    rt = dt/2;                 % throat radius
    re = de/2;                 % exit radius
    
    % Choose Rao-style defaults
    theta_n = deg2rad(theta_n);     % match (nozzle wall) angle
    theta_e = deg2rad(theta_e);      % exit angle
    Rc      = 1.5*rt;          % throat radius of curvature
    
    % Step 1: match point (end of circular arc)
    xm = Rc * sin(theta_n);
    rm = rt + Rc*(1 - cos(theta_n));
    b  = tan(theta_n);
    
    % Step 2: parabola length and coefficients
    Lp = (re - rm) / (0.5*tan(theta_e) + 0.5*tan(theta_n));
    a  = (tan(theta_e) - tan(theta_n)) / (2*Lp);
    
    % Total length
    L  = xm + Lp;
    
    % Build r(x)
    x  = linspace(0, L, points).';
    r  = zeros(size(x));
    % arc
    idx1 = x <= xm;
    r(idx1) = rt + Rc - sqrt(Rc^2 - x(idx1).^2);
    % parabola
    s = x(~idx1) - xm;
    r(~idx1) = a.*s.^2 + b.*s + rm;
end

    % figure()
    % plot(x,r)
