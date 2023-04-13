function [f, tau, state] = SoftConstraintsController(t, p, v, pose, omega, pd, psid, attFlag, enumR, enumA)
    global g m J JInv e3 I  enuml enumr dt;
    persistent Kp Kv Kg Komega Gammaf Gammatau varrhof varrhotau evhat para...
        kp kvarrtheta cp cvartheta Tp Tvartheta pcinfty pcu0 pcl0 varthetainfty vartheta0 cp0 ctheta0...
        sknd spnd psuf psu1df psu2df pslf psl1df psl2df svarthetand varthetasf varthetas1df varthetas2df Constr pcuf pclf varthetacf;
    if isempty(Kp)
        Kp = 5;
        Kv = 5;
        Kg = 20; 
        Komega = 20;
        Gammaf = 300;
        Gammatau = 300;
        varrhof = zeros(3,1);
        varrhotau = zeros(3,1);
        para = [500,500,0.1];
        
        
        kp = 4;
        kvarrtheta = 2;
        cp = 0.6;
        Tp = 2;
        Tvartheta = 1;
        pcinfty = 0.03;
        pcu0 = [0.3; 0.3; 0.3];
        pcl0 = [0.3; 0.3; 0.3];
        varthetainfty = 0.1;
        vartheta0 = 0.3;
        
        sknd =@(k,c0,T,tt,x0,xu0,xl0,n) 1 * double(n == 0) + double(tt<T & (x0>=xu0 | x0<=-xl0)) .* prod(k+2-n:k+1) ./ (-T)^n .* (c0 - 1) * (1 - tt/T)^(k+1-n);
        Constr =@(tt, n, init, infinity) infinity * double(n == 0) + (-1)^n * exp(-tt) * (init - infinity);
        pcuf =@(tt, n) Constr(tt, n, pcu0, pcinfty);
        pclf =@(tt, n) Constr(tt, n, pcl0, pcinfty);
        varthetacf =@(tt, n) Constr(tt, n, vartheta0, varthetainfty);
    end
    R = ToSO3(pose, attFlag);
    
    ptilde = p - pd(:, 1);
    if isempty(spnd)
        cp0 = ptilde ./ cp ./ min(pcu0, pcl0);
        spnd =@(tt,n) sknd(kp,cp0,Tp,tt,ptilde,pcu0,pcl0,n);
        psuf =@(tt) spnd(tt, 0) .* pcuf(tt, 0);
        psu1df =@(tt) spnd(tt, 1) .* pcuf(tt, 0) + spnd(tt, 0) .* pcuf(tt, 1);
        psu2df =@(tt) spnd(tt, 2) .* pcuf(tt, 0) + 2 * spnd(tt, 1) .* pcuf(tt, 1) + spnd(tt, 0) .* pcuf(tt, 2);
        pslf =@(tt) spnd(tt, 0) .* pclf(tt, 0);
        psl1df =@(tt) spnd(tt, 1) .* pclf(tt, 0) + spnd(tt, 0) .* pclf(tt, 1);
        psl2df =@(tt) spnd(tt, 2) .* pclf(tt, 0) + 2 * spnd(tt, 1) .* pclf(tt, 1) + spnd(tt, 0) .* pclf(tt, 2);
    end
    psl = pslf(t);
    psl1d = psl1df(t);
    psl2d = psl2df(t);
    psu = psuf(t);
    psu1d = psu1df(t);
    psu2d = psu2df(t);
    pcu = pcuf(t, 0);
    pcl = pclf(t, 0);
    
    ep = atanh((2 * ptilde - (psu - psl)) ./ (psu + psl));
    
    Dp = 0.5 * diag((psl + psu) .* (1 - tanh(ep).^2));
    DpInv = 2 * diag(1 ./ (psl + psu) ./ (1 - tanh(ep).^2));
    pdu = 0.5 * (psu1d .* (tanh(ep) + 1));
    pdl = 0.5 * (psl1d .* (tanh(ep) - 1));
    
    vtilde = v - pd(:, 2);
    vtilded = -Dp * Kp * ep + pdu + pdl;
    ev = vtilde - vtilded;
    ep1d = - Kp * ep + DpInv * ev;
    Dp1d = 0.5 * diag((psl + psu) .* (-2 * tanh(ep) .* (1 - tanh(ep).^2)) .* ep1d + (psl1d + psu1d) .* (1 - tanh(ep).^2));
    pdu1d = 0.5 * (psu1d .* (1 - tanh(ep).^2) .* ep1d + psu2d .* (tanh(ep) + 1));
    pdl1d = 0.5 * (psl1d .* (1 - tanh(ep).^2) .* ep1d + psl2d .* (tanh(ep) - 1));
    

    if isempty(evhat)
        evhat = zeros(3, 2);
        evhat(:,1) = ev;
    end
    ev1d = evhat(:, 2);
    for i = 1 : 3
        evhat(i,:) = Smoother(ev(i), evhat(i,:)', t, para)';
    end
    vtilded1d = -Dp1d * Kp * ep - Dp * Kp * ep1d + pdu1d + pdl1d;
    
    dfhat = varrhof + m * Gammaf * v;
    Fd = m * (Kv * ev + g * e3 - pd(:, 3) - vtilded1d + DpInv * ep) + dfhat;

    [gtilde, G, GInv, omegatilde, A, omegadbar, data] = GibbsVector(t, R, Fd, psid, omega, enumR, enumA, 'R');
    Fd = data{1};
    fd1d = data{2};
    fd = norm(Fd);
    f = Fd' * R(:, 3);
    varrhof = varrhof + dt * (ev / m - Gammaf * (dfhat + m * g * e3 - f * R(:, 3)));

    

    gtilde1d = G * omegatilde;
    gtilde2n = gtilde' * gtilde;
    
    varthetatilde = 2 * atan(norm(gtilde));

    if isempty(svarthetand)
        cvartheta = double(varthetatilde == 0) * 1 + double(varthetatilde ~= 0) * varthetatilde / 0.6;
        ctheta0 = varthetatilde/cvartheta/vartheta0;
        svarthetand =@(tt,n) sknd(kvarrtheta,ctheta0,Tvartheta,tt,varthetatilde,vartheta0,vartheta0,n);
        varthetasf =@(tt) svarthetand(tt, 0) * varthetacf(tt, 0);
        varthetas1df =@(tt) svarthetand(tt, 1) * varthetacf(tt, 0) + svarthetand(tt, 0) * varthetacf(tt, 1);
        varthetas2df =@(tt) svarthetand(tt, 2) * varthetacf(tt, 0) + 2 * svarthetand(tt, 1) * varthetacf(tt, 1) + svarthetand(tt, 0) * varthetacf(tt, 2);
    end
    varthetas = varthetasf(t); 
    varthetas1d = varthetas1df(t);
    varthetas2d = varthetas2df(t);
    varthetac = varthetacf(t, 0);
    cg = cos(varthetas/2)^2;
    cg1d = - 0.5 * varthetas1d * sin(varthetas);
    cg2d = - 0.5 * (varthetas2d * sin(varthetas) + varthetas1d^2 * cos(varthetas));
    if enumR == enuml
        Bf = I;
        Bf1d = zeros(3,3);
        b = R(:,3);
        b1d = R * R3so3(omega) * e3;
    elseif enumR == enumr
        Bf = R;
        Bf1d = R * R3so3(omega);
        b = e3;
        b1d = zeros(3,1);
    end
    eg = gtilde / (1 - cg * (1 + gtilde2n));
    Hp = 2 * fd * Bf * R3so3(b) * (gtilde' * b * R3so3(b) - I) / (1 + gtilde2n);
    Hp1d = 2 * (fd1d * Bf * R3so3(b) + fd * Bf1d * R3so3(b) + fd * Bf * R3so3(b1d)) * (gtilde' * b * R3so3(b) - I) / (1 + gtilde2n)...
           + 2 * fd * Bf * R3so3(b) * ((gtilde1d' * b * R3so3(b) + gtilde' * b1d * R3so3(b) + gtilde' * b * R3so3(b1d)) / (1 + gtilde2n) - (2 * gtilde' * gtilde1d * (gtilde' * b * R3so3(b) - I)) / (1 + gtilde2n)^2);
    H = (1 - cg * (1 + gtilde2n)) * Hp;

    H1d = (1 - cg * (1 + gtilde2n)) * Hp1d - (cg1d * (1 + gtilde2n) + 2 * cg * gtilde' * gtilde1d) * Hp;
    GInv1d = - (2 * (1 + gtilde2n) * enumA * R3so3(gtilde1d) + 4 * gtilde' * gtilde1d * (I - enumA * R3so3(gtilde))) / (1 + gtilde2n)^2;
    Omega = ((1 - cg * (1 + gtilde2n)) * I + 2 * cg * (gtilde * gtilde')) / (1 - cg * (1 + gtilde2n))^2 * G;
    eg1d = Omega * omegatilde + cg1d * (1 + gtilde2n) / (1 - cg * (1 + gtilde2n)) * eg;
    OmegaInv = GInv * (1 - cg * (1 + gtilde2n)) / (1 - cg * (1 - gtilde2n)) * ((1 - cg * (1 + gtilde2n)) * I - 2 * cg * R3so3(gtilde)^2);
    OmegaInv1d = (4 * cg * (cg - 1) * gtilde' * gtilde1d - 2 * cg1d * gtilde2n) / (1 - cg * (1 - gtilde2n))^2 * GInv * ((1 - cg * (1 + gtilde2n)) * I - 2 * cg * R3so3(gtilde)^2)...
                 +(1 - cg * (1 + gtilde2n)) / (1 - cg * (1 - gtilde2n)) * GInv1d * ((1 - cg * (1 + gtilde2n)) * I - 2 * cg * R3so3(gtilde)^2)...
                 +(1 - cg * (1 + gtilde2n)) / (1 - cg * (1 - gtilde2n)) * GInv * ((-cg1d * (1 + gtilde2n) - 2 * cg * gtilde' * gtilde1d) * I - 2 * cg1d * R3so3(gtilde)^2 - 2 * cg * (R3so3(gtilde) * R3so3(gtilde1d) + R3so3(gtilde1d) * R3so3(gtilde)));
             
    omegatilded = -OmegaInv * (Kg * eg + cg1d * (1 + gtilde2n)/(1-cg*(1+gtilde2n))*eg - H' * ev / m);
    eomega = omegatilde - omegatilded;

    
    omegatilded1d = -OmegaInv1d * (Kg * eg + cg1d * (1 + gtilde2n)/(1-cg*(1+gtilde2n)) * eg - H' * ev / m)...
                    -OmegaInv * (Kg * eg1d + cg1d * (1 + gtilde2n)/(1-cg*(1+gtilde2n)) * eg1d - H1d' * ev / m - H * ev1d / m ...
                    +((cg1d^2-cg*cg2d) * (1+gtilde2n)^2 + cg2d * (1+gtilde2n) + 2 * cg1d * gtilde' * gtilde1d) / (1-cg*(1+gtilde2n))^2 * eg);
    
    dtauhat = varrhotau + Gammatau * J * omega;
    tau = cross(omega, J*omega) - dtauhat + J * A' * (-Komega * eomega - omegadbar + omegatilded1d - Omega' * eg);
    varrhotau = varrhotau + dt * (JInv' * A' * eomega - Gammatau * (dtauhat - cross(omega,J*omega) + tau));

    state{1} = ptilde;
    state{2} = varthetatilde;
    state{3} = dfhat;
    state{4} = dtauhat;
    state{5} = pcu;
    state{6} = pcl;
    state{7} = psu;
    state{8} = psl;
    state{9} = varthetac;
    state{10} = varthetas;
    if t < 0
        f = m * g;
        tau = zeros(3,1);
    end
end
