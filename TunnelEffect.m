%% Parameter initialisation
repertoire = './';
executable = 'Exercice8';
input = 'configuration.in';
omega = 0.003;
% xL = -200;
% n = 14;
% L = 400;
% E0 = 2 * pi^2 * n^2 / L^2;
E0 = 0.0257;
ratioVE = 0:0.1:2;
param = 'delta';
paramval = sqrt(2 * ratioVE * E0)/omega;
nsimul = length(paramval);

%% Simulations 
for i = 1:nsimul
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
for i = 1:1
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
%     V = data(:,2);
    V0 = 0.5 * omega^2 * paramval(i)^2;
    xbar = delta - sqrt(2 * E(1)) / omega;
    
    figure('Name',[param,' = ',num2str(paramval(i)),'; ',...
        'p_0 = ',num2str(ratioVE(i)),'; ',...
        'p_1 = ',num2str(V0/E(1)),'; ',...
        'E_1 = ',num2str(E(1)),'; ',...
        'V_0 = ',num2str(V0),'; ',...
        'width',num2str(2 * xbar)])
    subplot(1,2,1)
    hold on
    pcolor(x,t,psi2)
    if V0 > E(1)
        plot([xbar,xbar],[min(t),max(t)],'w--')
        plot([-xbar,-xbar],[min(t),max(t)],'w--')
    end
    hold off
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
figure
xlimits = [paramval(1),paramval(end)];
ylimits = [min(pnshift),max(pnshift)];
plot(paramval,pnshift,'b.',(ylimits - p(2))/p(1),ylimits,'k',xlimits,[0.5 0.5],'r')
nhalf = (0.5 - p(2))/p(1);
title(['n_{min} = ',num2str(nhalf)])

%% Best n
cmd = sprintf('wsl %s%s %s %s=%.15g %s=%.15g %s=%.15g %s=%s', repertoire, executable, input,...
                                             'n', nhalf, 'delta', delta,...
                                             'x0', -delta,...
                                             'output', ['n=' num2str(nhalf)]);
disp(cmd); system(cmd);
%%
fichier = 'n=11.0641';
Analyse(fichier,delta,omega)

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
for i = [1 10 20 30]
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
                                             'n',3,...
                                             'output', 'TestBox');
    disp(cmd); system(cmd);
%%
Analyse('TestBox',0,omega)

%% Custom periodic
cmd = sprintf('wsl %s%s %s %s=%s %s=%s %s=%.15g %s=%s', repertoire, executable, input,...
                                             'potential','periodic',...
                                             'wavetype','eigenstate',...
                                             'n',3,...
                                             'output', 'TestPer');
    disp(cmd); system(cmd);
%%
Analyse('TestPer',0,omega)

%% Energy spectrum

