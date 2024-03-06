function [gtilde, G, GInv, omegatilde, A, omegadbar, data] = GibbsVector(t, pose, Fd, psid, omega, enumR, enumA, attFlag)
    global enuml enumr I;
    

    [desPose, Rd1d, Rd2d, data] = desirePose(t, Fd, psid, attFlag);
    R = ToSO3(pose, attFlag);
    Rd = ToSO3(desPose, attFlag);
    omegad = so3R3(Rd' * Rd1d);
    omegad1d = so3R3(Rd1d'*Rd1d + Rd'*Rd2d);
Rd1d'*Rd1d + Rd'*Rd2d
    if enumR == enuml
        Rtilde = R * Rd';
        if enumA == enuml
            omegatilde = R * (omega - omegad);
            A = R;
            omegadbar = -R * (omegad1d + cross(omega, omegad));
        elseif enumA == enumr
            omegatilde = Rd * (omega - omegad);
            A = Rd;
            omegadbar = -Rd * (omegad1d + cross(omega, omegad));
        end
    elseif enumR == enumr
        Rtilde = Rd' * R;
        if enumA == enuml
            omegatilde = Rtilde * omega - omegad;
            A = Rtilde;
            omegadbar = cross(omegatilde, omegad) - omegad1d;
        elseif enumA == enumr
            omegatilde = omega - Rtilde' * omegad;
            A = I;
            omegadbar = cross(omegatilde, omega) - Rtilde' * omegad1d;
        end
    end
	gtilde = so3R3(Rtilde - Rtilde') / (trace(Rtilde) + 1);
    G = 0.5 * (eye(3) + enumA * R3so3(gtilde) + gtilde * gtilde');
    GInv = 2 * (eye(3) - enumA * R3so3(gtilde)) / (1 + gtilde' * gtilde);
end

function [pose, Rd1d, Rd2d, data] = desirePose(t, Fd0, psid, attFlag)
    global dt;
    persistent Fd0hat sig QLast kq dQLast;
    if isempty(Fd0hat)
        Fd0hat = zeros(3,3);
        Fd0hat(:, 1) = Fd0;% 记录初值
    end

    Fd = Fd0hat(:, 1);
    Fd1d = Fd0hat(:, 2);
    Fd2d = Fd0hat(:, 3);
    for i = 1 : 3
        Fd0hat(i,:) = Smoother(Fd0(i), Fd0hat(i,:)', t, [500,500,0.0])';
    end
    fd = norm(Fd);
    
    bd = [cos(psid(1)); sin(psid(1)); 0];
    bd1d = psid(2) * [-sin(psid(1)); cos(psid(1)); 0];
    bd2d = psid(2)^2 * [-cos(psid(1)); -sin(psid(1)); 0] + psid(3) * [-sin(psid(1)); cos(psid(1)); 0];
    b3d = Fd / norm(Fd);
    b3d1d = - R3so3(b3d)^2 * Fd1d / fd;
    b3d2d = -((2*b3d1d*b3d'+b3d*b3d1d')*Fd1d + R3so3(b3d)^2*Fd2d)/fd;
    b2d = cross(b3d, bd);
    b2d = b2d / norm(b2d);
    b2d1d = - (R3so3(b2d)^2*(cross(b3d1d,bd)+cross(b3d,bd1d)))/norm(cross(b3d,bd));
    b2d2d = -((2*b2d1d*b2d'+b2d*b2d1d')*(cross(b3d1d,bd)+cross(b3d,bd1d)) + R3so3(b2d)^2*(cross(b3d2d,bd)+2*cross(b3d1d,bd1d)+cross(b3d,bd2d)))/norm(cross(b3d,bd));
    b1d = cross(b2d, b3d);
    b1d1d = cross(b2d1d,b3d) + cross(b2d,b3d1d);
    b1d2d = cross(b2d2d,b3d) + 2*cross(b2d1d,b3d1d) + cross(b2d,b3d2d);
    fd1d = b3d'*Fd1d;
    data{1} = Fd;
    data{2} = fd1d;
    Rd = [b1d, b2d, b3d];
    Rd1d = [b1d1d, b2d1d, b3d1d];
    Rd2d = [b1d2d, b2d2d, b3d2d];
    if attFlag == 'Q'
        if isempty(sig)
           sig = 1; 
           QLast = [1; 0; 0; 0];
           kq = 0.5;
	   dQLast = [1; 0; 0; 0];
        end
        pose = 0.5 * [sqrt(abs(1 + Rd(1,1) + Rd(2,2) + Rd(3,3)));
                      sign(Rd(3,2) - Rd(2,3)) * sqrt(abs(1 + Rd(1,1) - Rd(2,2) - Rd(3,3)));
                      sign(Rd(1,3) - Rd(3,1)) * sqrt(abs(1 - Rd(1,1) + Rd(2,2) - Rd(3,3)));
                      sign(Rd(2,1) - Rd(1,2)) * sqrt(abs(1 - Rd(1,1) - Rd(2,2) + Rd(3,3)))];
                  % 绝对值是为了防止由计算误差引起的微小负数带来的影响
        pose = sig * pose / norm(pose);
        
        % 变号检测  仅在特技飞行时可能用到，即当qd0可能为0时有必要
	dQ = (pose - QLast) / dt
 	ddQ = abs(dQ - dQLast)
        if min(diag(dQ) * dQLast) < 0 && max(ddQ) > kq
            sig = -sig;
	    pose = -pose;
     	    dQ = (pose - QLast) / dt
        end
        QLast = pose;
	dQLast = dQ;
        %结束
        
    elseif attFlag == 'R'
        pose = Rd;
    end
 
end
