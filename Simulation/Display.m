function Display
clear all;
close all;
clc;

global color t FontSize LineWide FigFontSize FigureCount;

color=[235 237 233;
       238 121 089;
       083 081 100;
       153 188 172;
       050 120 138;
       234 229 227;
       194 081 096;
       225 210 121;
       102 061 116;
       212 221 225;
       235 225 169;
       221 107 123;
       021 029 041;
       210 057 024;
       229 168 075;
       093 163 157;
       050 120 138;
       093 163 157;
       254 220 064;
       245 243 242;
       146 129 135;
       087 100 112;
       158 078 086;
       236 217 199]/255;

LineWide=1.8;
FontSize=11;
FigFontSize=11;

% for i = 1:24
%     plot(0:1, i*[1,1], 'color', color(i,:), 'linewidth', 11);
%     hold on;
%     axis('tight');
% end

d0=load('E.mat');
d1=load('PID.mat');
d2=load('LQ.mat');
d3=load('ADRC.mat');
d4=load('Soft.mat');
t=d0.t;

%% Dual
FigureCount = 10;
disp(['事件触发时长占比：',num2str(trapz(t,d0.varsigma)/d0.Ts)]);
disp(['事件触发次数：',num2str(d0.TriggerTimes)]);
myfigure();
for i = 1:3
    subplot(3,1,i);
    myplot([d0.pd(i,:);d0.pdhat(i,:);d0.p(i,:)], [5 14 13], {'-',':','--'}, ['p_',num2str(i),'/\textrm{m}']);
end
legend({'$$\textit{\textbf{p}}_d$$','$$\hat{\textit{\textbf{p}}}_d$$','$$\textit{\textbf{p}}$$'},'interpreter','latex','FontSize', FigFontSize);

myfigure();
subplot(3,1,1);
myplot(d0.Theta(1,:), 13, {'-'}, '\phi/\textrm{rad}');
subplot(3,1,2);
myplot(d0.Theta(2,:), 13, {'-'}, '\theta/\textrm{rad}');
subplot(3,1,3);
myplot([d0.psid(1,:);d0.psidhat(1,:);d0.Theta(3,:)], [5 14 13], {'-',':','--'}, '\psi/\textrm{rad}');
legend({'$$\psi_d$$','$$\hat{\psi}_d$$','$$\psi$$'},'interpreter','latex','FontSize', FigFontSize);

myfigure();
for i = 1:3
    subplot(3,1,i);
    myplot([d0.df(i,:);d0.dfhat(i,:)], [14 13], {'-','--'}, ['d_{f',num2str(i),'}/\textrm{N}']);
end
legend({'$$\textit{\textbf{d}}_f$$','$$\hat{\textit{\textbf{d}}}_{f}$$'},'interpreter','latex','FontSize', FigFontSize);

myfigure();
for i = 1:3
    subplot(3,1,i);
    myplot([d0.dtau(i,:);d0.dtauhat(i,:)], [14 13], {'-','--'}, ['d_{\tau',num2str(i),'}/\textrm{N}']);
end
legend({'$$\textit{\textbf{d}}_\tau$$','$$\hat{\textit{\textbf{d}}}_{\tau}$$'},'interpreter','latex','FontSize', FigFontSize);

myfigure();
for i = 1:4
    subplot(4,1,i);
    myplot(d0.T(i,:), 13, {'-'}, ['T_',num2str(i),'/\textrm{N}']);
end

myfigure();
subplot(3,1,1);
myplot(d0.e(1,:), 13, {'-'}, 'e_\textrm{total}');
subplot(3,1,2);
myplot(d0.varsigma(1,:), 13, {'-'}, '\varsigma');
subplot(3,1,3);
plot(t, d0.TriggeredPoint(1,:), 'x', 'color', color(13,:),'linewidth',0.5);
xlabel('$t/\textrm{s}$','interpreter','latex', 'FontSize', FontSize),ylabel('$(t_{k+1}-t_k)/\textrm{s}$','interpreter','latex', 'FontSize', 11);
axis('tight');

myfigure();
myplot3({d0.pd,d0.pdhat,d0.p},[5 14 13], {'-',':','--'});
legend({'$$\textit{\textbf{p}}_d$$','$$\hat{\textit{\textbf{p}}}_d$$','$$\textit{\textbf{p}}$$'},'interpreter','latex','FontSize', FigFontSize);

%% PID
FigureCount = 20;
myfigure();
for i = 1:3
    subplot(4,1,i);
    myplot(d1.ptilde(i,:), 13, {'-'}, ['\tilde{p}_',num2str(i),'/\textrm{m}']);
end
subplot(4,1,4);
myplot(d1.varthetatilde(1,:), 13, {'-'}, '|\tilde{\vartheta}|/\textrm{rad}');

myfigure();
subplot(4,1,1);
myplot(d1.f(1,:), 13, {'-'}, 'f/\textrm{N}');
for i = 1:3
    subplot(4,1,i+1);
    myplot(d1.tau(i,:), 13, {'-'}, ['\tau_',num2str(i),'/\textrm{N}\cdot\textrm{m}']);
end

myfigure();
for i = 1:3
    subplot(3,1,i);
    myplot([d1.pd(i,:);d1.pdhat(i,:);d1.p(i,:)], [5 14 13], {'-', ':', '--'}, ['p_',num2str(i),'/\textrm{m}']);
end
legend({'$$\textit{\textbf{p}}_d$$','$$\hat{\textit{\textbf{p}}}_d$$','$$\textit{\textbf{p}}$$'},'interpreter','latex','FontSize', FigFontSize);

myfigure();
myplot3({d1.pd,d1.pdhat,d1.p},[5 14 13], {'-','-','--'});
legend({'$$\textit{\textbf{p}}_d$$','$$\hat{\textit{\textbf{p}}}_d$$','$$\textit{\textbf{p}}$$'},'interpreter','latex','FontSize', FigFontSize);

%% LQ
FigureCount = 30;

myfigure();
for i = 1:3
    subplot(4,1,i);
    myplot(d2.ptilde(i,:), 13, {'-'}, ['\tilde{p}_',num2str(i),'/\textrm{m}']);
end
subplot(4,1,4);
myplot(d2.varthetatilde(1,:), 13, {'-'}, '|\tilde{\vartheta}|/\textrm{rad}');

myfigure();
subplot(4,1,1);
myplot(d2.f(1,:), 13, {'-'}, 'f/\textrm{N}');
for i = 1:3
    subplot(4,1,i+1);
    myplot(d2.tau(i,:), 13, {'-'}, ['\tau_',num2str(i),'/\textrm{N}\cdot\textrm{m}']);
end

myfigure();
for i = 1:3
    subplot(3,1,i);
    myplot([d2.pd(i,:);d2.pdhat(i,:);d2.p(i,:)], [5 14 13], {'-', ':', '--'}, ['p_',num2str(i),'/\textrm{m}']);
end
legend({'$$\textit{\textbf{p}}_d$$','$$\hat{\textit{\textbf{p}}}_d$$','$$\textit{\textbf{p}}$$'},'interpreter','latex','FontSize', FigFontSize);

myfigure();
myplot3({d2.pd,d2.pdhat,d2.p},[5 14 13], {'-','-','--'});
legend({'$$\textit{\textbf{p}}_d$$','$$\hat{\textit{\textbf{p}}}_d$$','$$\textit{\textbf{p}}$$'},'interpreter','latex','FontSize', FigFontSize);

%% ADRC
FigureCount = 40;
myfigure();
for i = 1:3
    subplot(4,1,i);
    myplot(d3.ptilde(i,:), 13, {'-'}, ['\tilde{p}_',num2str(i),'/\textrm{m}']);
end
subplot(4,1,4);
myplot(d3.varthetatilde(1,:), 13, {'-'}, '|\tilde{\vartheta}|/\textrm{rad}');

myfigure();
subplot(4,1,1);
myplot(d3.f(1,:), 13, {'-'}, 'f/\textrm{N}');
for i = 1:3
    subplot(4,1,i+1);
    myplot(d3.tau(i,:), 13, {'-'}, ['\tau_',num2str(i),'/\textrm{N}\cdot\textrm{m}']);
end

myfigure();
for i = 1:3
    subplot(3,1,i);
    myplot([d3.pd(i,:);d3.pdhat(i,:);d3.p(i,:)], [5 14 13], {'-', ':', '--'}, ['p_',num2str(i),'/\textrm{m}']);
end
legend({'$$\textit{\textbf{p}}_d$$','$$\hat{\textit{\textbf{p}}}_d$$','$$\textit{\textbf{p}}$$'},'interpreter','latex','FontSize', FigFontSize);

myfigure();
myplot3({d3.pd,d3.pdhat,d3.p},[5 14 13], {'-','-','--'});
legend({'$$\textit{\textbf{p}}_d$$','$$\hat{\textit{\textbf{p}}}_d$$','$$\textit{\textbf{p}}$$'},'interpreter','latex','FontSize', FigFontSize);

buf = FontSize;
buf1 = LineWide;
FontSize = 8;
LineWide = 1;
myfigure();
for i = 1:3
    for j = 1:3
        subplot(3,3,3*(i-1)+j)
        myplot([d3.zp(3*(i-1)+j,:);d3.zphat(3*(i-1)+j, :)], [14 13], {'-','--'}, ['z_{p',num2str(i),num2str(j),'}']);
    end
end
legend({'$$\textit{\textbf{z}}_{pi},i=1,2,3$$','$$\hat{\textit{\textbf{z}}}_{pi},i=1,2,3$$'},'interpreter','latex','FontSize', FigFontSize);

myfigure();
for i = 1:3
    for j = 1:3
        subplot(3,3,3*(i-1)+j)
        myplot([d3.zg(3*(i-1)+j,:);d3.zghat(3*(i-1)+j, :)], [14 13], {'-','--'}, ['z_{g',num2str(i),num2str(j),'}']);
    end
end
legend({'$$\textit{\textbf{z}}_{gi}, i=1,2,3$$','$$\hat{\textit{\textbf{z}}}_{gi}, i=1,2,3$$'},'interpreter','latex','FontSize', FigFontSize);
FontSize = buf;
LineWide = buf1;
%% Soft
FigureCount = 50;
myfigure();
for i = 1:3
    subplot(3,1,i);
    myplot([d4.df(i,:);d4.dfhat(i,:)], [14 13], {'-','--'}, ['d_{f',num2str(i),'}/\textrm{N}']);
end
legend({'$$\textit{\textbf{d}}_f$$','$$\hat{\textit{\textbf{d}}}_{f}$$'},'interpreter','latex','FontSize', FigFontSize);

myfigure();
for i = 1:3
    subplot(3,1,i);
    myplot([d4.dtau(i,:);d4.dtauhat(i,:)], [14 13], {'-','--'}, ['d_{\tau',num2str(i),'}/\textrm{N}']);
end
legend({'$$\textit{\textbf{d}}_\tau$$','$$\hat{\textit{\textbf{d}}}_{\tau}$$'},'interpreter','latex','FontSize', FigFontSize);

myfigure();
subplot(4,1,1);
myplot(d4.f(1,:), 13, {'-'}, 'f/\textrm{N}');
for i = 1:3
    subplot(4,1,i+1);
    myplot(d4.tau(i,:), 13, {'-'}, ['\tau_',num2str(i),'/\textrm{N}\cdot\textrm{m}']);
end


myfigure();
for i = 1:3
    subplot(3,1,i);
    myplot([d4.ptilde(i,:); d4.pcu(i, :); d4.psu(i,:); -d4.pcl(i,:); -d4.psl(i,:)], [13 5 14 5 14], {'-', '-', '--', '-', '--'}, ['\tilde{p}_',num2str(i),'/\textrm{m}']);
end
legend({'$$\tilde{\textit{\textbf{p}}}$$','$$\textit{\textbf{p}}_{u}(-\textit{\textbf{p}}_{l})$$','$$\textit{\textbf{p}}_{su}(-\textit{\textbf{p}}_{sl})$$'},'interpreter','latex','FontSize', FigFontSize);
myfigure();
myplot([d4.varthetatilde(1,:); d4.varthetac(1,:); d4.varthetacs(1,:)], [13 5 14], {'-','-','--'}, '|\tilde{\vartheta}|/\textrm{rad}');
legend({'$$|\tilde{\vartheta}|$$','$$\vartheta_{c}$$','$$\vartheta_s$$'},'interpreter','latex','FontSize', FigFontSize);

myfigure();
myplot3({d4.pd,d4.p},[14 13], {'-','--'});
legend({'$$\textit{\textbf{p}}_d$$','$$\textit{\textbf{p}}$$'},'interpreter','latex','FontSize', FigFontSize);

end



function myplot(v, coloridx, linetype, name)
    global color t FontSize LineWide;
    n = size(v, 1);
    for i = 1 : n
        plot(t, v(i, :), linetype{i}, 'color', color(coloridx(i), :),'linewidth',LineWide);
        hold on;
    end
    xlabel('$t/\textrm{s}$','interpreter','latex', 'FontSize', FontSize),ylabel(['$',name,'$'],'interpreter','latex', 'FontSize', FontSize);
    axis('tight');
end

function myplot3(v, coloridx, linetype)
    global color FontSize LineWide;
    n = size(v, 2);
    for i = 1 : n
        plot3(v{i}(1,:),v{i}(2,:),v{i}(3,:), linetype{i}, 'color', color(coloridx(i), :),'linewidth',LineWide);
        hold on;
    end
    xlabel('$x/\textrm{m}$','interpreter','latex', 'FontSize', FontSize),ylabel('$y/\textrm{m}$','interpreter','latex', 'FontSize', FontSize);zlabel('$z/\textrm{m}$','interpreter','latex', 'FontSize', FontSize);
    axis('tight');
end

function myfigure
global FigureCount;

FigureCount=FigureCount+1;
figure(FigureCount);
end

