clc;
close all;
clear all;


%% 调试参数
global debugFlag
debugFlag = false;
FigureCount=0;
%% 运行参数
global dt
testtime = 0.0;
Ts = 20;% s
times = (Ts + testtime) * 1000;
dt = Ts / times; %0.001s
length = times + 1;
t = linspace(-testtime, Ts, length);

DisturbanceSignalEnable = true;
%% 姿态及姿态误差选择

global enuml enumr;
enuml = 1;
enumr = -1;
attFlag = 'R';
enumR = enumr;
enumA = enumr;
%% 控制器选择
PIDFlag = 1;
LQFlag = 2;
ADRCFlag = 3;
SoftConstraintFlag = 4;

ControllerSelectFlag = SoftConstraintFlag;

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


if attFlag == 'R'
    R = zeros(3, 3, length);
    R(:, :, 1) = eye(3);
elseif attFlag == 'Q'
    Q = zeros(4, length);
    Q(1, 1) = 1;
end
pd = zeros(3, 5, length);
psid = zeros(5, length);
ptilde = zeros(3, length);
varthetatilde = zeros(1, length);
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

% % 轨迹  
for i =  1 : length
    pd(:,1, i)=[5*log(t(i)+1);
                cos(2.3*t(i));
                -sin(2.3*t(i))-2*log(t(i)+1)+0.];
    pd(:,2,i)=[5/(t(i)+1);
               -2.3*sin(2.3*t(i));
               -2.3*cos(2.3*t(i))-2/(t(i)+1)];
    pd(:,3,i)=[-5/(t(i)+1)^2;
               -2.3^2*cos(2.3*t(i));
               2.3^2*sin(2.3*t(i))+2/(t(i)+1)^2];
    pd(:,4,i)=[10/(t(i)+1)^3;
               2.3^3*sin(2.3*t(i));
               2.3^3*cos(2.3*t(i))-4/(t(i)+1)^3];
    psid(1,i) = -pi * cos(1.2*t(i));
    psid(2,i) = 1.2*pi * sin(1.2*t(i));
    psid(3,i) = 1.2^2*pi * cos(1.2*t(i));
    psid(4,i) = -1.2^3*pi * sin(1.2*t(i));
end
% 轨迹平滑
pdhat = zeros(3, 5, length);
psidhat = zeros(5, length);
para = [300,2,1];
    for i = 1 : times
        for n = 1 : 3
            pdhat(n, :, i+1) = Smoother(pd(n, 1, i),  pdhat(n, :, i)', t(i), para)';
        end
            psidhat(:, i+1) = Smoother(psid(1, i), psidhat(:, i), t(i), para);
    end

if ControllerSelectFlag == PIDFlag
    Controller = @PIDControrller;
elseif ControllerSelectFlag == LQFlag
    Controller = @LQControrller;
elseif ControllerSelectFlag == ADRCFlag
    Controller = @ADRController;
    zphat = zeros(9,length);
    zghat = zeros(9,length);
    zp = zeros(9, length);
    zg = zeros(9, length);
    zptilde = zeros(9, length);
    zgtilde = zeros(9, length);
elseif ControllerSelectFlag == SoftConstraintFlag
    Controller = @SoftConstraintsController;
    dfhat = zeros(3, length);
    dtauhat = zeros(3, length);
    pcu = zeros(3, length);
    pcl = zeros(3, length);
    psu = zeros(3, length);
    psl = zeros(3, length);
    varthetac = zeros(1, length);
    varthetacs = zeros(1, length);
    pd0 = [-0.2; 0.2; -0.8];
    pdhat(:,1,:) = pdhat(:,1,:) + pd0;
    psidhat(1, :) = psidhat(1, :);
    pd = pdhat;
    psid = psidhat;
end


for i = 1 : times

    if attFlag == 'R'
        [f(1,i), tau(:,i), state] = Controller(t(i), p(:, i), v(:, i), R(:,:,i), omega(:,i), pdhat(:,:,i), psidhat(:,i), attFlag, enumA, enumR);
        [p(:, i+1), v(:, i+1), R(:,:,i+1), omega(:,i+1)] = QuadrotorRotationMatrixModel(p(:,i), v(:,i), R(:,:,i), omega(:,i), f(1,i), tau(:,i), df(:,i), dtau(:,i));
    elseif attFlag == 'Q'
        [f(1,i), tau(:,i), state] = Controller(t(i), p(:, i), v(:, i), Q(:,i), omega(:,i), pdhat(:,:,i), psidhat(:,i), attFlag, enumA, enumR);
        [p(:, i+1), v(:, i+1), Q(:,i+1), omega(:,i+1)] = QuadrotorQuaternionModel(p(:,i), v(:,i), Q(:,i), omega(:,i), f(:,i), tau(:,i), df(:,i), dtau(:,i));
    end
    ptilde(:,i) = state{1};
    varthetatilde(:,i) = state{2};
    if ControllerSelectFlag == ADRCFlag
        zphat(:,i+1) = state{3};
        zghat(:,i+1) = state{4};
        zp(:, i) = [state{1};
                    v(:, i) - pdhat(:, 2, i);
                    df(:, i)/m - pdhat(:,3,i)];
        zg(:, i) = [state{5};
                    state{6};
                    state{8}*(-cross(omega(:,i),J*omega(:,i))+dtau(:,i))+state{7}+2*state{5}'*state{6}*state{6}/(1+state{5}'*state{5})];        
    elseif ControllerSelectFlag == SoftConstraintFlag
        dfhat(:, i) = state{3};
        dtauhat(:, i) = state{4};
        pcu(:, i) = state{5};
        pcl(:, i) = state{6};
        psu(:, i) = state{7};
        psl(:, i) = state{8};
        varthetac(:, i) = state{9};
        varthetacs(:, i) = state{10};
    end
end
i = length;
if attFlag == 'R'
    [f(1,i), tau(:,i), state] = Controller(t(i), p(:, i), v(:, i), R(:,:,i), omega(:,i), pdhat(:,:,i), psidhat(:,i), attFlag, enumA, enumR);
elseif attFlag == 'Q'
    [f(1,i), tau(:,i), state] = Controller(t(i), p(:, i), v(:, i), Q(:,i), omega(:,i), pdhat(:,:,i), psidhat(:,i), attFlag, enumA, enumR);
end
ptilde(:,i) = state{1};
varthetatilde(:,i) = state{2};
if ControllerSelectFlag == ADRCFlag
    zp(:, i) = [state{1};
                v(:, i) - pdhat(:, 2, i);
                df(:, i)/m - pdhat(:,3,i)];
    zg(:, i) = [state{5};
                state{6};
                state{8}*(-cross(omega(:,i),J*omega(:,i))+dtau(:,i))+state{7}+2*state{5}'*state{6}*state{6}/(1+state{5}'*state{5})];    
    if debugFlag
    FigureCount=FigureCount+1;
    figure(FigureCount);
    subplot(3,3,1);
        plot(t, zp(1, :), t, zphat(1, :),'r--');
    subplot(3,3,2);
        plot(t, zp(2, :), t, zphat(2, :),'r--');
    subplot(3,3,3);
        plot(t, zp(3, :), t, zphat(3, :),'r--');
    subplot(3,3,4);
        plot(t, zp(4, :), t, zphat(4, :),'r--');
    subplot(3,3,5);
        plot(t, zp(5, :), t, zphat(5, :),'r--');
    subplot(3,3,6);
        plot(t, zp(6, :), t, zphat(6, :),'r--');
    subplot(3,3,7);
        plot(t, zp(7, :), t, zphat(7, :),'r--');
    subplot(3,3,8);
        plot(t, zp(8, :), t, zphat(8, :),'r--');
    subplot(3,3,9);
        plot(t, zp(9, :), t, zphat(9, :),'r--');

     FigureCount=FigureCount+1;
    figure(FigureCount);
    subplot(3,3,1);
        plot(t, zg(1, :), t, zghat(1, :),'r--');
    subplot(3,3,2);
        plot(t, zg(2, :), t, zghat(2, :),'r--');
    subplot(3,3,3);
        plot(t, zg(3, :), t, zghat(3, :),'r--');
    subplot(3,3,4);
        plot(t, zg(4, :), t, zghat(4, :),'r--');
    subplot(3,3,5);
        plot(t, zg(5, :), t, zghat(5, :),'r--');
    subplot(3,3,6);
        plot(t, zg(6, :), t, zghat(6, :),'r--');
    subplot(3,3,7);
        plot(t, zg(7, :), t, zghat(7, :),'r--');
    subplot(3,3,8);
        plot(t, zg(8, :), t, zghat(8, :),'r--');
    subplot(3,3,9);
        plot(t, zg(9, :), t, zghat(9, :),'r--');
    end
elseif ControllerSelectFlag == SoftConstraintFlag
    dfhat(:, i) = state{3};
    dtauhat(:, i) = state{4};
    pcu(:, i) = state{5};
    pcl(:, i) = state{6};
    psu(:, i) = state{7};
    psl(:, i) = state{8};
    varthetac(:, i) = state{9};
    varthetacs(:, i) = state{10};
    
    
end
pdhat = reshape(pdhat(:,1,:),size(p));
pd = reshape(pd(:, 1, :), size(p));

if ControllerSelectFlag == PIDFlag
    save('PID');
elseif ControllerSelectFlag == LQFlag
    save('LQ');
elseif ControllerSelectFlag == ADRCFlag
    save('ADRC');
elseif ControllerSelectFlag == SoftConstraintFlag
    save('Soft');
end


if debugFlag

FontSize = 8;
FigFontSize = 10;
LineWide = 1;

if ControllerSelectFlag == SoftConstraintFlag
FigureCount=FigureCount+1;
figure(FigureCount);
subplot(4,1,1);
plot(t,ptilde(1,:),'k', t, pcu(1, :), 'b', t, psu(1, :), 'r--', t, -pcl(1, :), 'b', t, -psl(1, :), 'r--');axis('tight');
xlabel('$t$','interpreter','latex', 'FontSize', 10),ylabel('$\tilde{p}_{1}$','interpreter','latex', 'FontSize', 10);
subplot(4,1,2);
plot(t,ptilde(2,:),'k', t, pcu(2, :), 'b', t, psu(2, :), 'r--', t, -pcl(2, :), 'b', t, -psl(2, :), 'r--');axis('tight');
xlabel('$t$','interpreter','latex', 'FontSize', 10),ylabel('$\tilde{p}_{2}$','interpreter','latex', 'FontSize', 10);
subplot(4,1,3);
plot(t,ptilde(3,:),'k', t, pcu(3, :), 'b', t, psu(3, :), 'r--', t, -pcl(3, :), 'b', t, -psl(3, :), 'r--');axis('tight');
xlabel('$t$','interpreter','latex', 'FontSize', 10),ylabel('$\tilde{p}_{3}$','interpreter','latex', 'FontSize', 10);
subplot(4,1,4);
plot(t,varthetatilde(1,:),'k', t, varthetac(1, :), 'b', t, varthetacs(1, :), 'r--');axis('tight');
xlabel('$t$','interpreter','latex', 'FontSize', 10),ylabel('$|\tilde{\vartheta}|$','interpreter','latex', 'FontSize', 10);
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
end

FigureCount=FigureCount+1;
figure(FigureCount);
subplot(3,1,1);
    plot(t, pd(1, :),'r--', t, p(1, :));
subplot(3,1,2);
    plot(t, pd(2, :),'r--', t, p(2, :));
subplot(3,1,3);
    plot(t, pd(3, :),'r--', t, p(3, :));
    
FigureCount=FigureCount+1;
figure(FigureCount);
plot3(pd(1,:),pd(2,:),pd(3,:),'g-.',p(1,:),p(2,:),p(3,:),'b','linewidth',LineWide);
legend({'$$p_d$$','$$p$$'},'interpreter','latex','FontSize', FigFontSize);
xlabel('$x$ [m]','interpreter','latex', 'FontSize', FontSize),ylabel('$y$ [m]','interpreter','latex', 'FontSize', FontSize),zlabel('$z$ [m]','interpreter','latex', 'FontSize', FontSize);


FigureCount=FigureCount+1;
figure(FigureCount);
plot(t,f);
xlabel('$t$','interpreter','latex', 'FontSize', 10),ylabel('$f$','interpreter','latex', 'FontSize', 10);

FigureCount=FigureCount+1;
figure(FigureCount);
subplot(3,1,1);
plot(t,tau(1,:));
xlabel('$t$','interpreter','latex', 'FontSize', 10),ylabel('$\tau_1$','interpreter','latex', 'FontSize', 10);
subplot(3,1,2);
plot(t,tau(2,:));
xlabel('$t$','interpreter','latex', 'FontSize', 10),ylabel('$\tau_2$','interpreter','latex', 'FontSize', 10);
subplot(3,1,3);
plot(t,tau(3,:));
xlabel('$t$','interpreter','latex', 'FontSize', 10),ylabel('$\tau_3$','interpreter','latex', 'FontSize', 10);

FigureCount=FigureCount+1;
figure(FigureCount);
subplot(4,1,1);
plot(t,ptilde(1,:));axis('tight');
xlabel('$t$','interpreter','latex', 'FontSize', 10),ylabel('$\tilde{p}_{1}$','interpreter','latex', 'FontSize', 10);
subplot(4,1,2);
plot(t,ptilde(2,:));axis('tight');
xlabel('$t$','interpreter','latex', 'FontSize', 10),ylabel('$\tilde{p}_{2}$','interpreter','latex', 'FontSize', 10);
subplot(4,1,3);
plot(t,ptilde(3,:));axis('tight');
xlabel('$t$','interpreter','latex', 'FontSize', 10),ylabel('$\tilde{p}_{3}$','interpreter','latex', 'FontSize', 10);
subplot(4,1,4);
plot(t,varthetatilde(1,:));axis('tight');
xlabel('$t$','interpreter','latex', 'FontSize', 10),ylabel('$|\tilde{\vartheta}|$','interpreter','latex', 'FontSize', 10);


% FigureCount=FigureCount+1;
% figure(FigureCount);
% subplot(3,1,1);
% plot(t,df(1,:));
% xlabel('$t$','interpreter','latex', 'FontSize', 10),ylabel('$d_{f1}$','interpreter','latex', 'FontSize', 10);
% subplot(3,1,2);
% plot(t,df(2,:));
% xlabel('$t$','interpreter','latex', 'FontSize', 10),ylabel('$d_{f2}$','interpreter','latex', 'FontSize', 10);
% subplot(3,1,3);
% plot(t,df(3,:));
% xlabel('$t$','interpreter','latex', 'FontSize', 10),ylabel('$d_{f3}$','interpreter','latex', 'FontSize', 10);

end


