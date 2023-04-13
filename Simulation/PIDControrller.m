function [f, tau, state] = PIDControrller(t, p, v, pose, omega, pd, psid, attFlag, enumR, enumA)
    global dt g e3 m J;
    persistent Kpp Kpi Kpd Kgp Kgi Kgd ptildeInt gtildeInt;
    if isempty(gtildeInt)
        omegac = 20;
        omegacg = 10;
        Kpp = omegac^2;
        Kpi = 0;
        Kpd = 2 * omegac;
        Kgp = omegacg^2;
        Kgi = 0;
        Kgd = 2 * omegacg;
        ptildeInt = zeros(3,1);
        gtildeInt = zeros(3,1);
    end
    ptilde = p - pd(:,1);
    vtilde = v - pd(:,2);
    ptildeInt = ptildeInt + dt * ptilde;
    FdLump = -Kpp * ptilde - Kpi * ptildeInt - Kpd * vtilde;
    Fd = -m * (FdLump - g * e3 + pd(:,3)); 
%     Fd = -(norm(ptilde)+norm(vtilde)+1)*(FdLump - g * e3 + pd(:,3));
%     f = Fd' * (ToSO3(pose, attFlag) * e3);
    
    [gtilde, G, GInv, omegatilde, A, omegadbar, data] = GibbsVector(t, pose, Fd, psid, omega, enumR, enumA, attFlag);
    Fd = data{1};
    f = Fd' * (ToSO3(pose, attFlag) * e3);
    
    gtilde1d = G * omegatilde;
    gtildeInt = gtildeInt + gtilde * dt;
    tauLump = -Kgp * gtilde - Kgi * gtildeInt - Kgd * gtilde1d;
    tau = J * A' * (GInv * (tauLump - 2 * gtilde' * gtilde1d * gtilde1d / (1 + gtilde' * gtilde)) - omegadbar) + cross(omega, J * omega);
%     tau = 0.1 * A' * (GInv * (tauLump - 2 * gtilde' * gtilde1d * gtilde1d / (1 + gtilde' * gtilde)) - omegadbar);

	varthetatilde = 2 * atan(norm(gtilde));
    state{1} = ptilde;
    state{2} = varthetatilde;
end