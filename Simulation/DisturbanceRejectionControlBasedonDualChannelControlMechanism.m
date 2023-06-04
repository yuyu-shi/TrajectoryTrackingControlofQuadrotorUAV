
function DisturbanceRejectionControlBasedonDualChannelControlMechanism

clear all;

clear all;
clear all;
clear all;

close all;
clc;
%% 调试参数
global debugFlag
debugFlag = true;
FigureCount=0;
%% 运行参数
global dt
Ts = 20;% s
times = 20000;
dt = Ts / times; %0.001s
length = times + 1;
t = linspace(0, Ts, length);

DisturbanceSignalEnable = true;
%% 飞行器参数
global m J JInv g e3 C CInv cf Te I;
cf = 1;
d = 0.315;
ctau = 8.004e-4;
m = 4.34;
I = eye(3);
g = 9.81;
J = diag([0.0820, 0.0845, 0.1377]);
JInv = J^-1;
e3 = [0, 0, 1]';
C = [cf cf cf cf;0 -d*cf 0 d*cf; d*cf 0 -d*cf 0; ctau -ctau ctau -ctau];
CInv=[1/4/cf 0 1/2/d/cf 1/4/ctau;1/4/cf -1/2/d/cf 0 -1/4/ctau; 1/4/cf 0 -1/2/d/cf 1/4/ctau; 1/4/cf 1/2/d/cf 0 -1/4/ctau];
Te = m * g / 4 / cf;

%% 状态量
p = zeros(3, length);
v = zeros(3, length);
omega = zeros(3, length);
df = zeros(3, length);
dtau = zeros(3, length);
f = zeros(1, length);
tau = zeros(3, length);
T = zeros(4, length);
Theta = zeros(3, length);

pd = zeros(3, 5, length);
psid = zeros(1, 5, length);


if DisturbanceSignalEnable
    for i = 1 : length
        df(:,i)=[(cos(10*t(i))-sin(9*t(i))-cos(8*t(i))-sin(7*t(i)));
                 (cos(20*t(i))-sin(8*t(i))-cos(5*t(i))-sin(9*t(i)));
                 (-cos(11*t(i))-sin(5*t(i))+cos(6*t(i))-sin(8*t(i)))];
        dtau(:,i)=0.1*[(sin(11*t(i))+cos(12*t(i))+sin(13*t(i))-cos(7*t(i)));
                   (sin(11*t(i))-cos(3*t(i))+sin(13*t(i))+cos(10*t(i)));
                   (sin(12*t(i))+cos(12*t(i))+sin(7*t(i))-cos(7*t(i)))];
    end
    
end

%% 控制器
% 双通道控制
%参数定义
    Tbreak = 7000;
%变量
    pdhat = zeros(3, 5, length);
    psidhat = zeros(1, 5, length);
    dfhat = zeros(3, length);
    dtauhat = zeros(3, length);
    TriggeredPoint = zeros(1, length) * NaN;
    varsigma = zeros(1, length);
    e = zeros(1, length);
% 轨迹  
    for ti = 1 : length
        pd(:,1,ti) = [3.5 / pi * cos((t(ti)) * 0.2 * pi);
                      3.5 / pi * sin((t(ti)) * 0.2 * pi);
                      -1.5 * log(t(ti) + 1) - 1];
        pd(:,2,ti) = [-0.7 * sin((t(ti)) * 0.2 * pi);
                      0.7 * cos((t(ti)) * 0.2 * pi);
                      -1.5 * (t(ti) + 1)^-1];
        pd(:,3,ti) = [-0.14 * pi * cos((t(ti)) * 0.2 * pi);
                      -0.14 * pi * sin((t(ti)) * 0.2 * pi);
                      1.5 * (t(ti) + 1)^-2]; 
        pd(:,4,ti) = [0.028 * pi^2 * sin((t(ti)) * 0.2 * pi);
                      -0.028 * pi^2 * cos((t(ti)) * 0.2 * pi);
                      -3 * (t(ti) + 1)^-3];
        pd(:,5,ti) = [0.0056 * pi^3 * cos((t(ti)) * 0.2 * pi);
                      0.0056 * pi^3 * sin((t(ti) - t(Tbreak)) * 0.2 * pi);
                      9 * (t(ti) - t(Tbreak) + 1)^-4];
        psid(1,1,ti) = pi * sin(t(ti)) / 10 - 1;
        psid(1,2,ti) = pi * cos(t(ti)) / 10;
        psid(1,3,ti) = -pi * sin(t(ti)) / 10;
        psid(1,4,ti) = -pi * cos(t(ti)) / 10;
        psid(1,5,ti) = pi * sin(t(ti)) / 10;
    end
% 轨迹平滑
para = [300,3.5,1];
% para = [300,10,5];
    for i = 1 : times
        for n = 1 : 3
            pdhat(n, :, i+1) = Smoother(pd(n, 1, i),  pdhat(n, :, i)', (i - 1) * dt, para)';
        end
        psidhat(:, :, i+1) = Smoother(psid(:, 1, i), psidhat(:, :, i)', (i - 1) * dt, para)';
    end

% 轨迹跟踪               
    for i = 1 : times 
        [f(:, i), tau(:, i), varsigma(1, i), dfhat(:, i), dtauhat(:, i), TriggeredPoint(:, i)] = DualChannelController(p(:, i), v(:, i), Theta(:, i), omega(:, i), pdhat(:, :, i), psidhat(:, :, i));
        T(:, i) = CInv * [f(:, i); tau(:, i)];
        [p(:, i+1), v(:, i+1), Theta(:, i+1), omega(:, i+1)] = QuadrotorEulerAngleModel(p(:, i), v(:, i), Theta(:, i), omega(:, i), f(:, i), tau(:, i), df(:, i), dtau(:, i));
    end
    i = length;
    [f(:, i), tau(:, i), varsigma(1, i), dfhat(:, i), dtauhat(:, i), TriggeredPoint(:, i)] = DualChannelController(p(:, i), v(:, i), Theta(:, i), omega(:, i), pdhat(:, :, i), psidhat(:, :, i));
    T(:, i) = CInv * [f(:, i); tau(:, i)];
    
    pd = reshape(pd(:,1,:),size(p));
    psid = reshape(psid(:,1,:),size(Theta(3,:)));
    pdhat = reshape(pdhat(:,1,:),size(p));
    psidhat = reshape(psidhat(:,1,:),size(Theta(3,:)));
    for i = 1 : length
    	e(i) = norm([p(:,i)-pdhat(:,i); Theta(3,i)-psidhat(1,i)]);
    end
    TriggeredPoint(1) = 0;
    TriggerTimes = 1;
    last = 0;
    for i = 2 : length
        if(varsigma(i-1)) == 0
                last = t(i);
        end
        if varsigma(i) == 1 && TriggeredPoint(i) == 1
            TriggerTimes = TriggerTimes + 1;
            TriggeredPoint(i) = t(i) - last;
            last = t(i);
        end
    end     

save('E');

if debugFlag

FontSize = 8;
FigFontSize = 10;
LineWide = 1;
FigureCount=FigureCount+1;
figure(FigureCount);
plot(t, e, t, 0.016*varsigma);

FigureCount=FigureCount+1;
figure(FigureCount);
plot(t,TriggeredPoint,'x',t, 0.02*varsigma);
% hold on;
% for j=1:length
%   if ~isnan(TriggeredPoint(j))
%     line([t(j),t(j)],[0,TriggeredPoint(j)]); hold on;
%   end
% end

FigureCount=FigureCount+1;
figure(FigureCount);
subplot(3,1,1);
plot(t,df(1,:),'b',t,dfhat(1,:),'r--','linewidth',LineWide);
xlabel('$t$','interpreter','latex', 'FontSize', FontSize),ylabel('$d_{f1}$','interpreter','latex', 'FontSize', FontSize);
subplot(3,1,2);
plot(t,df(2,:),'b',t,dfhat(2,:),'r--','linewidth',LineWide);
xlabel('$t$','interpreter','latex', 'FontSize', FontSize),ylabel('$d_{f2}$','interpreter','latex', 'FontSize', FontSize);
subplot(3,1,3);
plot(t,df(3,:),'b',t,dfhat(3,:),'r--','linewidth',LineWide);
xlabel('$t$','interpreter','latex', 'FontSize', FontSize),ylabel('$d_{f3}$','interpreter','latex', 'FontSize', FontSize);
legend({'$$d_f$$','$$\hat{d}_f$$'},'interpreter','latex','FontSize', FigFontSize);

FigureCount=FigureCount+1;
figure(FigureCount);
subplot(3,1,1);
plot(t,dtau(1,:),'b',t,dtauhat(1,:),'r--','linewidth',LineWide);
xlabel('$t$','interpreter','latex', 'FontSize', FontSize),ylabel('$d_{\tau 1}$','interpreter','latex', 'FontSize', FontSize);
subplot(3,1,2);
plot(t,dtau(2,:),'b',t,dtauhat(2,:),'r--','linewidth',LineWide);
xlabel('$t$','interpreter','latex', 'FontSize', FontSize),ylabel('$d_{\tau 2}$','interpreter','latex', 'FontSize', FontSize);
subplot(3,1,3);
plot(t,dtau(3,:),'b',t,dtauhat(3,:),'r--','linewidth',LineWide);
xlabel('$t$','interpreter','latex', 'FontSize', FontSize),ylabel('$d_{\tau 3}$','interpreter','latex', 'FontSize', FontSize);
legend({'$$d_\tau$$','$$\hat{d}_\tau$$'},'interpreter','latex','FontSize', FigFontSize);





FigureCount=FigureCount+1;
figure(FigureCount);
subplot(3,1,1);
    plot(t, pd(1, :),'r--', t, pdhat(1, :), t, p(1, :));
subplot(3,1,2);
    plot(t, pd(2, :),'r--', t, pdhat(2, :), t, p(2, :));
subplot(3,1,3);
    plot(t, pd(3, :),'r--', t, pdhat(3, :), t, p(3, :));
FigureCount=FigureCount+1;
figure(FigureCount);
subplot(3,1,1);
    plot(t, Theta(1, :));
subplot(3,1,2);
    plot(t, Theta(2, :));
subplot(3,1,3);
    plot(t, Theta(3, :), t, psid(1, :), t, psidhat(1, :));  

    
    
FigureCount=FigureCount+1;
figure(FigureCount);
plot3(pd(1,:),pd(2,:),pd(3,:),'r--',pdhat(1,:),pdhat(2,:),pdhat(3,:),'g-.',p(1,:),p(2,:),p(3,:),'b','linewidth',LineWide);
legend({'$$p_d$$','$$\hat{p}_d$$','$$p$$'},'interpreter','latex','FontSize', FigFontSize);
xlabel('$x$ [m]','interpreter','latex', 'FontSize', FontSize),ylabel('$y$ [m]','interpreter','latex', 'FontSize', FontSize),zlabel('$z$ [m]','interpreter','latex', 'FontSize', FontSize);
FigureCount=FigureCount+1;
figure(FigureCount);
plot(t,T(1,:),'k');



FigureCount=FigureCount+1;
figure(FigureCount);
subplot(4,1,1);
plot(t,T(1,:),'k');
axis('tight');
xlabel('$t$ [s]','interpreter','latex', 'FontSize', FigFontSize),ylabel('$f_{1}[N]$','interpreter','latex', 'FontSize', FigFontSize);
subplot(4,1,2);
plot(t,T(2,:),'k');
axis('tight');
xlabel('$t$ [s]','interpreter','latex', 'FontSize', FigFontSize),ylabel('$f_{2}[N]$','interpreter','latex', 'FontSize', FigFontSize);
subplot(4,1,3);
plot(t,T(3,:),'k');
axis('tight');
xlabel('$t$ [s]','interpreter','latex', 'FontSize', FigFontSize),ylabel('$f_{3}[N]$','interpreter','latex', 'FontSize', FigFontSize);
subplot(4,1,4);
plot(t,T(4,:),'k');
axis('tight');
xlabel('$t$ [s]','interpreter','latex', 'FontSize', FigFontSize),ylabel('$f_{4}[N]$','interpreter','latex', 'FontSize', FigFontSize);

% FigureCount=FigureCount+1;
% figure(FigureCount);
% plot(t,f);
% xlabel('$t$','interpreter','latex', 'FontSize', 10),ylabel('$f$','interpreter','latex', 'FontSize', 10);

% FigureCount=FigureCount+1;
% figure(FigureCount);
% subplot(3,1,1);
% plot(t,tau(1,:));
% xlabel('$t$','interpreter','latex', 'FontSize', 10),ylabel('$F_1$','interpreter','latex', 'FontSize', 10);
% subplot(3,1,2);
% plot(t,tau(2,:));
% xlabel('$t$','interpreter','latex', 'FontSize', 10),ylabel('$F_2$','interpreter','latex', 'FontSize', 10);
% subplot(3,1,3);
% plot(t,tau(3,:));
% xlabel('$t$','interpreter','latex', 'FontSize', 10),ylabel('$F_3$','interpreter','latex', 'FontSize', 10);

% FigureCount=FigureCount+1;
% figure(FigureCount);
% plot(t, e);
% FigureCount=FigureCount+1;
% figure(FigureCount);
% plot(t, varsigma);
% FigureCount=FigureCount+1;
% figure(FigureCount);
% plot(t, TriggerFlag);
% FigureCount=FigureCount+1;
% figure(FigureCount);
% plot(t, p(1, :), t, p(2, :), t, p(3, :));
% FigureCount=FigureCount+1;
% figure(FigureCount);
% subplot(3,1,1);
% plot(t,Theta(1,:),'k');
% xlabel('$t$','interpreter','latex', 'FontSize', 10),ylabel('$\phi$','interpreter','latex', 'FontSize', 10);
% subplot(3,1,2);
% plot(t,Theta(2,:),'k');
% xlabel('$t$','interpreter','latex', 'FontSize', 10),ylabel('$\theta$','interpreter','latex', 'FontSize', 10);
% subplot(3,1,3);
% plot(t,Theta(3,:),'k');
% xlabel('$t$','interpreter','latex', 'FontSize', 10),ylabel('$\psi$','interpreter','latex', 'FontSize', 10);
end
end

%% 基于双通道控制机制的粗轨迹跟踪控制
function [f, tau, varsigma, dfhat, dtauhat, TriggerFlag] = DualChannelController(p, v, Theta, omega, pd, psid)
    global he0 Omegae deltae0 C CInv J JInv g m e3 dt I Te;
    persistent epsilon Gammatau Gammaf  Kp Kv Kl Komega k1 k2 k3 S32 S23 I32 varrhotau varrhof T ethresholdup ethresholddown sigma;
    if isempty(epsilon) 
        ethresholdup = 0.012;
        ethresholddown = 0.01;
        he0 = 5;
        Omegae = 30;
        deltae0 = 0.2;
        epsilon = 4;
        Gammatau = 50;
        Gammaf = 50;
        Kp = 8;
        Kv = 8;
        Kl = 8;
        Komega = 8;
        k1 = 1;
        k2 = 1; 
        k3 = 0.1;
        S32 = [1 0 0;
               0 1 0];
        S23 = S32';
        I32 = diag([1,1,0]);
        varrhotau = [0 0 0]';
        varrhof = [0 0 0]';
        T = [Te Te Te Te]';
        sigma = 0;
    end
% 控制器参数
    phi = Theta(1);
    theta = Theta(2);
    psi = Theta(3);
	b3 = [cos(psi) * sin(theta) * cos(phi) + sin(psi) * sin(phi);
         sin(psi) * sin(theta) * cos(phi) - cos(psi) * sin(phi);
         cos(phi) * cos(theta)];
    R33 = b3(3);
    W = [1, tan(theta) * sin(phi), tan(theta) * cos(phi);
         0, cos(phi), -sin(phi);
         0, sin(phi) * sec(theta), cos(phi) * sec(theta)];
    E = R33^-1 * [sin(psi) * sec(phi), sin(psi) * tan(theta) * tan(phi) + cos(psi) * sec(theta), 0;
                  -cos(psi) * sec(phi), -cos(psi) * tan(theta) * tan(phi) + sin(psi) * sec(theta), 0;
                  0, sin(phi) * cos(phi), cos(phi)^2];
    EInv = R33 * [-cos(psi) * sin(theta) * sin(phi) + sin(psi) * cos(phi), -sin(psi) * sin(theta) * sin(phi) - cos(psi) * cos(phi), 0;
                  cos(psi) * cos(theta), sin(psi) * cos(theta), 0;
                  -cos(psi) * tan(phi) * cos(theta), -sin(psi) * tan(phi) * cos(theta), sec(phi)^2];
    e11=-sin(psi) * sin(phi) - cos(psi) * sin(theta) * cos(phi);
    e12=-cos(psi) * cos(theta) * sin(phi);
    e13=cos(psi) * cos(phi) + sin(psi) * sin(theta) * sin(phi);
    e14=cos(psi) * sin(phi) - sin(psi) * sin(theta) * cos(phi);
    e15=-sin(psi) * cos(theta) * sin(phi);
    e16=sin(psi) * cos(phi) - cos(psi) * sin(theta) * sin(phi);
    e22=-cos(psi) * sin(theta);
    e23=-sin(psi) * cos(theta);
    e25=-sin(psi) * sin(theta);
    e26=cos(psi) * cos(theta);
    e31=-cos(psi) * cos(theta) * sec(phi)^2;
    e32=cos(psi) * sin(theta) * tan(phi);
    e33=sin(psi) * cos(theta) * tan(phi);
    e34=-sin(psi) * cos(theta) * sec(phi)^2;
    e35=sin(psi) * sin(theta) * tan(phi);
    e36=-cos(psi) * tan(phi) * cos(theta);
    e37=2*sec(phi)^2 * tan(phi);
    Theta1d = W * omega;
    EInv1d = R33 * [[e11 e12 e13;
                     0 e22 e23;
                     e31 e32 e33] * Theta1d, [e14 e15 e16;
                                            0 e25 e26;
                                            e34 e35 e36] * Theta1d, [0 0 0;
                                                                     0 0 0;
                                                                     e37 0 0] * Theta1d] + [-sin(phi)*cos(theta), -cos(phi)*sin(theta), 0] * Theta1d * R33^-1 * EInv;
                                                                 
    ptilde = p - pd(:,1);
    alpha1 = -Kp * ptilde + pd(:, 2);
    vtilde = v - alpha1;
    ptilde1d = v - pd(:, 2);
    alpha11d = -Kp * (v - pd(:,2)) + pd(:, 3);
    dfhat = varrhof + m * Gammaf * v;  %扰动观测器计算
    Fd = m * (ptilde + Kv * vtilde + g * e3 + dfhat / m - alpha11d);
    fd = Fd(3) / R33;
    alpha2 = S32 * Fd / Fd(3);
    l = [b3(1:2) / R33; psi];
    ld = [alpha2; psid(1,1)];
    ltilde = l - ld;
    
    
    Nf = -Fd(3)^-2 * S32 * R3so3(e3) * R3so3(Fd);
    Nf_f = Kv + Kp;
    Ndf = Nf_f + Gammaf;
    vtilde1df = -ptilde - Kv * vtilde - m^-1 * Fd(3) * I32 * ltilde;
    alpha12df = - Kp * (vtilde1df + alpha11d - pd(:,3)) + pd(:, 4);
    Fd1df = m * (v - pd(:, 2) + Kv * vtilde1df - alpha12df);
    alpha21df = - Fd(3)^-2 * S32 * R3so3(e3) * R3so3(Fd) * Fd1df;
    ld1df = [alpha21df; psid(1, 2)];
    alpha3f = EInv * (Fd(3) / m * I32 * vtilde - Kl * ltilde + ld1df - k1 * S23 * Nf * Ndf * Nf' * S32 * ltilde);
    
    omegatilde = omega - alpha3f;
   
    ltilde1df = Fd(3) / m * I32 * vtilde - Kl * ltilde + E * omegatilde;    
    vtilde2dff = -ptilde1d - Kv * vtilde1df - m^-1 * Fd1df(3) * I32 * ltilde - m^-1 * Fd(3) * I32 * ltilde1df;
    alpha13dff = -Kp * (vtilde2dff + alpha12df - pd(:, 4)) + pd(:, 5);
    Fd2dff = m * (vtilde1df + alpha11d - pd(:,3) + Kv * vtilde2dff - alpha13dff);
    Nfd1d = -Fd(3)^-2 * S32 * R3so3(e3) * R3so3(Fd1df) - 2 * e3' * Fd1df * Fd(3)^-1 * Nf;
    alpha22dff = Nf * Fd2dff - 2 * Fd1df(3) * Fd(3)^-1 * Nf * Fd1df;
    ld2dff = [alpha22dff; psid(1, 3)];
    alpha31dff = EInv1d * (Fd(3) / m * I32 * vtilde - Kl * ltilde + ld1df - k1 * S23 * Nf * Ndf * Nf' * S32 * ltilde)... 
                 + EInv * (Fd1df(3) / m * I32 * vtilde + Fd(3) / m * I32 * vtilde1df - Kl * ltilde1df + ld2dff...
                 - k1 * S23 * ((Nfd1d * Ndf * Nf' + Nf * Ndf * Nfd1d') * S32 * ltilde + Nf * Ndf * Nf' * S32 * ltilde1df));      
              
    Nff = -Nf * Nf_f * (I32 * ltilde * e3' - Fd(3) * S23 * Nf) - 2 * Fd(3)^-1 * Nf * Fd1df * e3' + Fd(3)^-2 * S32 * R3so3(e3) * R3so3(Fd1df);
    Nvff = m^-2 * Fd(3) *  EInv * I32 + EInv * S23 * Nf * (I - Kp * Kp - Nf_f * Kv);
    Nfff = EInv * (m^-1 * I32 * vtilde * e3' + Kl * S23 * Nf + S23 * Nff + k1 * S23 * Nf * Ndf * (Nf' * Nf)...
           + k1 * S23 * (Fd(3)^-2 * S32 * R3so3(e3) * R3so3(Ndf * Nf' * S32 *ltilde) - 2 * Fd(3)^-1 * Nf * Ndf * Nf' * S32 * ltilde * e3'));
    
    mutau = - k2 * Nfff * Ndf * Nfff' * omegatilde - k3 * (Nf_f * Nf_f') * omegatilde;
    
    dtauhat = varrhotau + Gammatau * J * omega; %扰动观测器计算
    
    taud = J * (-E' * ltilde - Komega * omegatilde + alpha31dff + mutau) + cross(omega, J * omega) - dtauhat;
    
    if norm([p-pd(:,1); Theta(3)-psid(1,1)]) < ethresholddown %精度足够
        sigma = 1;
    elseif norm([p-pd(:,1); Theta(3)-psid(1,1)]) > ethresholdup %精度不够
        sigma = 0;
    end
    varsigma = sigma;
    h = [(omegatilde' * (Nvff + Nfff * Nf_f) + ltilde' * S23 * Nf * Nf_f - vtilde' / m) * b3; JInv' * omegatilde];
    Td0 = CInv * [fd; taud];    
    
    Td = Td0 - varsigma * (deltae0 * Omegae + he0) * tanh(epsilon^-1 * (deltae0 * Omegae + he0) * C' * h);

    [T, TriggerFlag] = DualChannelMechanism(T, Td, varsigma);
    f = C(1, :) * T;
    tau = C(2:4,:) * T;
    varrhof = dt * ( - Gammaf * (dfhat + m * g * e3 - f * b3)) + varrhof; %扰动观测器更新
    varrhotau = dt * (JInv' * omegatilde - Gammatau * (dtauhat - cross(omega, J * omega) + tau)) + varrhotau; %扰动观测器更新
end
%带有事件触发的双通道机制
function [T, TriggerFlag] = DualChannelMechanism(TLast, Td, varsigma)
    global he0 Omegae deltae0 Te;
 
    TriggerFlag = NaN;
    if varsigma
        Tdelta = TLast - Te;
        if norm(Tdelta,'inf') < Omegae
            deltae = deltae0;
            he = he0;
        else
            deltae = 0;
            he = deltae0 * Omegae + he0;
        end
        if max(abs(TLast - Td) - deltae * abs(Tdelta)) < he
            T = TLast;
        else
            T = Td;
            TriggerFlag = 1;
        end
    else
        T = Td;
    end

end



