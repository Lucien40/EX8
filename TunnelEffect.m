%% Parameter initialisation
repertoire = './';
executable = 'Exercice8';
input = 'configuration.in';
omega = 0.003;
% xL = -200;
% n = 14;
% L = 400;
% E0 = 2 * pi^2 * n^2 / L^2;
E0 = 0.0258;
ratioVE = 0:0.1:2;
param = 'delta';
paramval = sqrt(2 * ratioVE * E0)/omega;
nsimul = length(paramval);

%% Simulations 
for i = [8,11,14]
    cmd = sprintf('wsl %s%s %s %s=%.15g %s=%.15g %s=%s', repertoire, executable, input,...
                                             param, paramval(i), 'x0', -paramval(i),...
                                             'output', [param '=' num2str(paramval(i))]);
    disp(cmd); system(cmd);
end

%% Wave analysis
numsim = 1;
fichier = [param '=' num2str(paramval(numsim))];
delta = paramval(numsim);
% psi2 = load([fichier '_psi2.out']);
% data = load([fichier '_obs.out']);
% t = data(:,1);
% probn = data(:,2);
% probp = data(:,3);
% E = data(:,4);
% data = load([fichier '_pot.out']);
% x = data(:,1);
% V = data(:,2);
disp('Data loaded')

Analyse(fichier,delta,omega)
% m = 100;
% liveWave(x,t,psi2,m)

%% All simulations
for i = [8,11,14]
    fichier = [param '=' num2str(paramval(i))];
    delta = paramval(i);
    psi2 = load([fichier '_psi2.out']);
    data = load([fichier '_obs.out']);
    t = data(:,1);
    probn = data(:,2);
    probp = data(:,3);
    E = data(:,4);
    data = load([fichier '_pot.out']);
    x = data(:,1);
    V = data(:,2);
    V0 = 0.5 * omega^2 * paramval(i)^2;
    xbar = delta - sqrt(2 * E(1)) / omega;
    
    H = 5;
    W = 8;
    figuA = figure;
    figuA.PaperUnits = 'centimeters';
    figuA.Units = 'centimeters';
    figuA.InvertHardcopy = 'on';
    figuA.PaperSize = [W H];
    figuA.PaperPosition = [0 0 W H];
    figuA.Position = [10 10 W H];
    
%     figure('Name',[param,' = ',num2str(paramval(i)),'; ',...
%         'p_0 = ',num2str(ratioVE(i)),'; ',...
%         'p_1 = ',num2str(V0/E(1)),'; ',...
%         'E_1 = ',num2str(E(1)),'; ',...
%         'V_0 = ',num2str(V0),'; ',...
%         'width',num2str(2 * xbar)])
%     subplot(1,2,1)
    hold on
    pcolor(x,t,psi2)
    if V0 >= E(1)
        plot([xbar,xbar],[min(t),max(t)],'w--','LineWidth',1)
        plot([-xbar,-xbar],[min(t),max(t)],'w--','LineWidth',1)
    end
    hold off
    shading interp
    colormap jet
    c = colorbar;
    xlabel('x [m]')
    ylabel('t [s]')
    ylabel(c,'|\psi|^2')
    title(['$\Delta$ = ',num2str(paramval(i)),';   ',...
        '$V_0$ = ',num2str(V0),';   ',...
        '$p$ = ',num2str(V0/E(1))])
    
%     subplot(1,2,2)
    H = 5;
    W = 8;
    figuA = figure;
    figuA.PaperUnits = 'centimeters';
    figuA.Units = 'centimeters';
    figuA.InvertHardcopy = 'on';
    figuA.PaperSize = [W H];
    figuA.PaperPosition = [0 0 W H];
    figuA.Position = [10 10 W H];
    hold on
    plot(t,probn,'Color',colors(6,:),'LineWidth',1)
    plot(t,probp,'Color',colors(2,:),'LineWidth',1)
    plot(t,probn + probp,'Color',colors(10,:),'LineWidth',1)
    hold off
    grid on
    leg = legend('$P(x<0)$','$P(x>0)$','$P_\mathrm{total}$');
%     leg.Color = 'none';
    leg.Location = 'best';
    title(['$\Delta$ = ',num2str(paramval(i)),';   ',...
        '$V_0$ = ',num2str(V0),';   ',...
        '$p$ = ',num2str(V0/E(1))])
    xlabel('t [s]')
    ylabel('probability')
end

%% Transmission probability 1/2
repertoire = './';
executable = 'Exercice8';
input = 'configuration.in';
omega = 0.003;
xL = -200;
delta = 64;
L = 400;
param = 'n';
paramval = 1:40;
% paramval = 10:0.1:12;
nsimul = length(paramval);
vg = 2 * pi * paramval / L;
% tshift = dist ./ vg;
% indshift = round(tshift + 1);

%% Simulations 
for i = 1:nsimul
    cmd = sprintf('wsl %s%s %s %s=%.15g %s=%.15g %s=%.15g %s=%s', repertoire, executable, input,...
                                             param, paramval(i), 'delta', delta,...
                                             'x0', -delta,...
                                             'output', [param '=' num2str(paramval(i))]);
    disp(cmd); system(cmd);
end

%%
difflim = 1e-5;
% plainSize = 100;
indshift = zeros(nsimul,1);
for i = 6:40
    fichier = [param '=' num2str(paramval(i))];
    psi2 = load([fichier '_psi2.out']);
    data = load([fichier '_obs.out']);
    t = data(:,1);
    probn = data(:,2);
    probp = data(:,3);
%     E = data(:,4);
    data = load([fichier '_pot.out']);
    x = data(:,1);
    
    figure('Name',[param,' = ',num2str(paramval(i))])
    subplot(1,2,1)
    pcolor(x,t,psi2)
    shading interp
    colormap jet
    c = colorbar;
    xlabel('x [m]')
    ylabel('t [s]')
    ylabel(c,'|\psi|^2')
    
    dp = diff(probn);
    index = (1:length(dp))';
    stable = abs(dp) < difflim;
    unstable = ~stable;
    plainNum = 0;
    cursor = 1;
    lb = 1;
    rb = 1;
    plainStop = 2;
%     if stable(1)
%         plainStop = 2;
%     else
%         plainStop = 1;
%     end
    while plainNum < plainStop
        ilim = index(index >= cursor & stable);
        lb = ilim(1);
        ilim = index(index >= cursor & unstable);
        rb = ilim(1) - 1;
        intSize = rb - lb;
        plainNum = plainNum + 1;
%         if intSize >= plainSize
%             plainNum = plainNum + 1;
%         end
        ilim = index(index >= rb + 1 & stable);
        cursor = ilim(1);
    end
    indshift(i) = round(mean([lb,rb]));
    
    subplot(1,2,2)
    hold on
    plot(t,probn,'r.',...
        t,probp,'b.',t,probn + probp,'m')
    plot([t(indshift(i)) t(indshift(i))],[0 1],'g')
    plot(t(1:end-1),stable,'k--')
    hold off
    xlabel('t [s]')
    ylabel('probability')
end
disp Finish

%%
pnshift = zeros(nsimul,1);
for i = 1:5
    fichier = [param '=' num2str(paramval(i))];
    data = load([fichier '_obs.out']);
    probn = data(:,2);
%     probp = data(:,3);
    pnshift(i) = probn(1) - probn(1500);
end
for i = 6:nsimul
    fichier = [param '=' num2str(paramval(i))];
    data = load([fichier '_obs.out']);
    probn = data(:,2);
%     probp = data(:,3);
    pnshift(i) = probn(1) - probn(indshift(i));
end
p = polyfit([11,12],[pnshift(11),pnshift(12)],1);
H = 5;
W = 8;
figuA = figure;
figuA.PaperUnits = 'centimeters';
figuA.Units = 'centimeters';
figuA.InvertHardcopy = 'on';
figuA.PaperSize = [W H];
figuA.PaperPosition = [0 0 W H];
figuA.Position = [10 10 W H];
xlimits = [paramval(1),paramval(end)];
ylimits = [min(pnshift),max(pnshift)];
hold on
plot(xlimits,[0.5 0.5],'Color',colors(6,:))
plot((ylimits - p(2))/p(1),ylimits,'Color',colors(1,:))
plot(paramval,pnshift,'.','Color',colors(2,:))
hold off
grid on
legend('$\Delta P = 0.5$','Linear regression','$\Delta P(n)$')
xlabel('n')
ylabel('$\Delta P$')
nhalf = (0.5 - p(2))/p(1);
title(['$n_{half}$ = ',num2str(nhalf)])

%% Best n
cmd = sprintf('wsl %s%s %s %s=%.15g %s=%.15g %s=%.15g %s=%s', repertoire, executable, input,...
                                             'n', nhalf, 'delta', delta,...
                                             'x0', -delta,...
                                             'output', ['n=' num2str(nhalf)]);
disp(cmd); system(cmd);
%%
fichier = 'n=11.0641';
% Analyse(fichier,delta,omega)

H = 5;
    W = 8;
    figuA = figure;
    figuA.PaperUnits = 'centimeters';
    figuA.Units = 'centimeters';
    figuA.InvertHardcopy = 'on';
    figuA.PaperSize = [W H];
    figuA.PaperPosition = [0 0 W H];
    figuA.Position = [10 10 W H];
    
%     figure('Name',[param,' = ',num2str(paramval(i)),'; ',...
%         'p_0 = ',num2str(ratioVE(i)),'; ',...
%         'p_1 = ',num2str(V0/E(1)),'; ',...
%         'E_1 = ',num2str(E(1)),'; ',...
%         'V_0 = ',num2str(V0),'; ',...
%         'width',num2str(2 * xbar)])
%     subplot(1,2,1)
    hold on
    pcolor(x,t,psi2)
%     if V0 >= E(1)
%         plot([xbar,xbar],[min(t),max(t)],'w--','LineWidth',1)
%         plot([-xbar,-xbar],[min(t),max(t)],'w--','LineWidth',1)
%     end
    hold off
    shading interp
    colormap jet
    c = colorbar;
    xlabel('x [m]')
    ylabel('t [s]')
    ylabel(c,'|\psi|^2')
%     title(['$\Delta$ = ',num2str(paramval(i)),';   ',...
%         '$V_0$ = ',num2str(V0),';   ',...
%         '$p$ = ',num2str(V0/E(1))])
    
%     subplot(1,2,2)
    H = 5;
    W = 8;
    figuA = figure;
    figuA.PaperUnits = 'centimeters';
    figuA.Units = 'centimeters';
    figuA.InvertHardcopy = 'on';
    figuA.PaperSize = [W H];
    figuA.PaperPosition = [0 0 W H];
    figuA.Position = [10 10 W H];
    hold on
    plot(t,probn,'Color',colors(6,:),'LineWidth',1)
    plot(t,probp,'Color',colors(2,:),'LineWidth',1)
    plot(t,probn + probp,'Color',colors(10,:),'LineWidth',1)
    hold off
    grid on
    leg = legend('$P(x<0)$','$P(x>0)$','$P_\mathrm{total}$');
%     leg.Color = 'none';
    leg.Location = 'best';
%     title(['$\Delta$ = ',num2str(paramval(i)),';   ',...
%         '$V_0$ = ',num2str(V0),';   ',...
%         '$p$ = ',num2str(V0/E(1))])
    xlabel('t [s]')
    ylabel('probability')

%% Box, periodic potential
repertoire = './';
executable = 'Exercice8';
input = 'configuration.in';
boxrange = 100;
omega = 0.003;
V0 = 0.1;
npot = 40;
param = 'n';
paramval = 1:30;
nsimul = length(paramval);
% potential type: Box, Per
pottype = 'Per';

%% Box Simulations 
for i = 1:nsimul
    cmd = sprintf('wsl %s%s %s %s=%s %s=%.15g %s=%.15g %s=%s', repertoire, executable, input,...
                                             'potential','box',...
                                             'boxrange',boxrange,...
                                             param, paramval(i),...
                                             'output', [param 'Box=' num2str(paramval(i))]);
    disp(cmd); system(cmd);
end

%% Periodic Simulations 
for i = 1:nsimul
    cmd = sprintf('wsl %s%s %s %s=%s %s=%.15g %s=%.15g %s=%.15g %s=%s', repertoire, executable, input,...
                                             'potential','periodic',...
                                             'V0',V0,...
                                             'npot',npot,...
                                             param, paramval(i),...
                                             'output', [param 'Per=' num2str(paramval(i))]);
    disp(cmd); system(cmd);
end

%% Wave analysis
numsim = 30;
fichier = [param pottype '=' num2str(paramval(numsim))];
delta = paramval(numsim);
psi2 = load([fichier '_psi2.out']);
data = load([fichier '_obs.out']);
t = data(:,1);
probn = data(:,2);
probp = data(:,3);
E = data(:,4);
data = load([fichier '_pot.out']);
x = data(:,1);
V = data(:,2);
disp('Data loaded')

Analyse(fichier,delta,omega)

%% All simulations
for i = 1:30 %[1 10 20 30]
    fichier = [param pottype '=' num2str(paramval(i))];
    delta = paramval(i);
    psi2 = load([fichier '_psi2.out']);
    data = load([fichier '_obs.out']);
    t = data(:,1);
    probn = data(:,2);
    probp = data(:,3);
    E = data(:,4);
    data = load([fichier '_pot.out']);
    x = data(:,1);
    V = data(:,2);
    
    figure('Name',[param,' = ',num2str(paramval(i)),'; ',...
        'E_1 = ',num2str(E(1))])
    subplot(1,2,1)
    pcolor(x,t,psi2)
    shading interp
    colormap jet
    c = colorbar;
    xlabel('x [m]')
    ylabel('t [s]')
    ylabel(c,'|\psi|^2')
    
    subplot(1,2,2)
    plot(t,probn,'r',...
        t,probp,'b',t,probn + probp,'m','LineWidth',1)
    xlabel('t [s]')
    ylabel('probability')
end


%% Custom box
cmd = sprintf('wsl %s%s %s %s=%s %s=%s %s=%.15g %s=%.15g %s=%.15g %s=%s', repertoire, executable, input,...
                                             'potential','box',...
                                             'wavetype','eigenstate',...
                                             'boxrange',100,...
                                             'infpot',0.001,...
                                             'n',1,...
                                             'output', 'TestBox');
    disp(cmd); system(cmd);
%%
boxr = 100;
dx = x(2) - x(1);
C = -0.5 / dx^2;
nx = length(x);
V0 = 0.001;
% V0 = 100000;
% Vbox = zeros(nx,1);
% Vbox(x <= -boxr | x >= boxr) = V0;
Vbox = 500 * sin(40 * 2.0 * pi * x / 400);
A = diag(-2*C*ones(nx,1) + Vbox,0) + diag(C*ones(nx-1,1),1) + diag(C*ones(nx-1,1),-1);
[Vec,D] = eig(A);
En = diag(D);
H = 5;
W = 8;
figuA = figure;
figuA.PaperUnits = 'centimeters';
figuA.Units = 'centimeters';
figuA.InvertHardcopy = 'on';
figuA.PaperSize = [W H];
figuA.PaperPosition = [0 0 W H];
figuA.Position = [10 10 W H];
hold on
for i = 1:5
    plot(x,Vec(:,i),'Color',colors(2*i,:))
end
% for i = 1:4
%     plot(x,Vec(:,149 + i),'Color',colors(2*i,:))
% end
% plot([boxr,boxr],[min(t),max(t)],'k--','LineWidth',1)
% plot([-boxr,-boxr],[min(t),max(t)],'k--','LineWidth',1)
xlabel('x [m]')
ylabel('$\Phi_n$')
legend('n = 1','n = 2','n = 3','n = 4','n = 5')
% legend('n = 150','n = 151','n = 152','n = 153')
hold off
grid on

H = 5;
W = 8;
figuA = figure;
figuA.PaperUnits = 'centimeters';
figuA.Units = 'centimeters';
figuA.InvertHardcopy = 'on';
figuA.PaperSize = [W H];
figuA.PaperPosition = [0 0 W H];
figuA.Position = [10 10 W H];
% subplot(2,1,1)
% hold on
% for i = 1:5
%     plot(i,En(i),'.','Color',colors(2*i,:))
% end
% xlabel('n')
% ylabel('$E_n$')
% hold off
% grid on
% box off
% subplot(2,1,2)
hold on
plot(En,'.','Color',colors(1,:))
xlabel('n')
ylabel('$E_n$')
hold off
grid on
box off

%%
for i = 1:301
    figure
    plot(x,Vec(:,i))
end

%%
fichier = 'TestBox';
    psi2 = load([fichier '_psi2.out']);
    data = load([fichier '_obs.out']);
    t = data(:,1);
    probn = data(:,2);
    probp = data(:,3);
%     E = data(:,4);
    data = load([fichier '_pot.out']);
    x = data(:,1);

H = 10;
    W = 8;
    figuA = figure;
    figuA.PaperUnits = 'centimeters';
    figuA.Units = 'centimeters';
    figuA.InvertHardcopy = 'on';
    figuA.PaperSize = [W H];
    figuA.PaperPosition = [0 0 W H];
    figuA.Position = [10 10 W H];
    
    hold on
    pcolor(x,t,psi2)
        plot([100,100],[min(t),max(t)],'w--','LineWidth',1)
        plot([-100,-100],[min(t),max(t)],'w--','LineWidth',1)
    hold off
    shading interp
    colormap jet
    c = colorbar;
    xlabel('x [m]')
    ylabel('t [s]')
    ylabel(c,'|\psi|^2')
    
    %%
Analyse('TestBox',0,omega)

%% Custom periodic
cmd = sprintf('wsl %s%s %s %s=%s %s=%.15g %s=%s', repertoire, executable, input,...
                                             'potential','periodic',...
                                             'n',20,...
                                             'output', 'TestPer');
    disp(cmd); system(cmd);
%%
fichier = 'nPer=30';
    psi2 = load([fichier '_psi2.out']);
    data = load([fichier '_obs.out']);
    t = data(:,1);
    probn = data(:,2);
    probp = data(:,3);
    meanx = data(:,5);
    meanp = data(:,7);
    meanp2 = data(:,8);
    data = load([fichier '_pot.out']);
    x = data(:,1);

H = 10;
    W = 8;
    figuA = figure;
    figuA.PaperUnits = 'centimeters';
    figuA.Units = 'centimeters';
    figuA.InvertHardcopy = 'on';
    figuA.PaperSize = [W H];
    figuA.PaperPosition = [0 0 W H];
    figuA.Position = [10 10 W H];
    
    hold on
    pcolor(x,t,psi2)
    hold off
    shading interp
    colormap jet
    c = colorbar;
    xlabel('x [m]')
    ylabel('t [s]')
    ylabel(c,'|\psi|^2')
    
    H = 5;
    W = 8;
    figuA = figure;
    figuA.PaperUnits = 'centimeters';
    figuA.Units = 'centimeters';
    figuA.InvertHardcopy = 'on';
    figuA.PaperSize = [W H];
    figuA.PaperPosition = [0 0 W H];
    figuA.Position = [10 10 W H];
    plot(t,meanx)
    xlabel('t [s]')
    ylabel('$\langle x\rangle$')
    
    H = 5;
    W = 8;
    figuA = figure;
    figuA.PaperUnits = 'centimeters';
    figuA.Units = 'centimeters';
    figuA.InvertHardcopy = 'on';
    figuA.PaperSize = [W H];
    figuA.PaperPosition = [0 0 W H];
    figuA.Position = [10 10 W H];
    plot(t,meanp)
    xlabel('t [s]')
    ylabel('$\langle p\rangle$')
    
%%
Analyse('TestPer',0,omega)

%% Energy spectrum

