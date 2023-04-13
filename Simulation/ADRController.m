function [f, tau, state] = ADRController(t, p, ~, pose, omega, pd, psid, attFlag, enumR, enumA)
    global g m J JInv e3;
    persistent kp1 kp2 kg1 kg2 ap1 ap2 ag1 ag2 zphat zghat zphat0 zghat0 uf0 utau0;
    if isempty(kp1)
        omegac = 20;
        kp1 = omegac^2;
        kp2 = 2 * omegac;
        kg1 = omegac^2;
        kg2 = 2 * omegac;
        ap1 = 1;
        ap2 = 1;
        ag1 = 1;
        ag2 = 1;
        zphat = zeros(9,1);
        zghat = zeros(9,1);
    end
    R = ToSO3(pose, attFlag);
    ptilde = p - pd(:,1);
    Fd0 = -kp1 * fal(zphat(1:3), ap1) - kp2 * fal(zphat(4:6), ap2);
    Fd = -m * (Fd0 - g * e3 - zphat(7:9));
    
    [gtilde, G, GInv, omegatilde, A, omegadbar, data] = GibbsVector(t, R, Fd, psid, omega, enumR, enumA, 'R');
    Fd = data{1};
    f = Fd' * R(:, 3);
    
    if isempty(zphat0)
        zphat0 = [ptilde; zeros(6, 1)];
    end
    zphat = ESO3order(ptilde - zphat0(1:3, :), -f*R(:,3)/m+g*e3, zphat) + zphat0;%观测下一时刻的扩张状态
    
    tau0 = -kg1 * fal(zghat(1:3), ag1) - kg2 * fal(zghat(4:6), ag2);
    tau = J * A' * GInv * (tau0 - zghat(7:9));
    
    if isempty(zghat0)
        zghat0 = [gtilde; zeros(6, 1)];
    end
    zghat = ESO3order(gtilde - zghat0(1:3,:), G*A*JInv*tau, zghat) + zghat0;%观测下一时刻的扩张状态
%     zghat = ESO3order(t, gtilde, G*A*JInv*tau, zghat);%观测下一时刻的扩张状态
    
    varthetatilde = 2 * atan(norm(gtilde));
    state{1} = ptilde;
    state{2} = varthetatilde;
    state{3} = zphat;
    state{4} = zghat;
    state{5} = gtilde;
    state{6} = G * omegatilde;
    state{7} = G * omegadbar;
    state{8} = G*A*JInv;
end
function zhat = ESO3order(z,B0u,zhat)
    global dt;
    persistent alpha01 alpha02 beta01 beta02 beta03;
    if isempty(alpha01)        
        alpha01 = 1;
        alpha02 = 1;
        omegao = 400;
        beta01 = 3 * omegao;
        beta02 = 3 * omegao^2;
        beta03 = omegao^3;
    end

    ztilde = z - zhat(1:3,1);
    zhat = zhat + dt * [zhat(4:6,1) + beta01 * ztilde;
                        zhat(7:9,1) + B0u + beta02 * fal(ztilde, alpha01);
                        beta03 * fal(ztilde, alpha02)];

end
function fe = fal(x, a)
    persistent delta;
    if isempty(delta)
        delta = 0.1;
    end
    c = double(x < delta);
    fe = c .* (x ./ delta^(1-a)) + (1 - c) .* (sign(x) .* abs(x).^a);
%     fe = x;
end

