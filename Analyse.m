function Analyse(fichier,delta,omega)
    %% Chargement des resultats %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    psi2 = load([fichier '_psi2.out']);
    data = load([fichier '_obs.out']);
    t = data(:,1);
    probn = data(:,2);
    probp = data(:,3);
    E = data(:,4);
    data = load([fichier '_pot.out']);
    x = data(:,1);

    %% Figures %%
    %%%%%%%%%%%%%
    figure('Name',['Analyse de ' fichier])
    subplot(2,2,1)
    plot(t,probn,'r',t,probp,'b',t,probn + probp,'m')
    grid
    xlabel('t [s]')
    ylabel('probability')

    subplot(2,2,2)
    V0 = 0.5 * omega^2 * delta^2;
    plot(t,E,'b')
    grid
    xlabel('t [s]')
    ylabel('E [J]')
    title(['E_1 = ',num2str(E(1)),'; V_0 = ',num2str(V0),...
        '; E_1 - V_0 = ',num2str(E(1)-V0)])

    subplot(2,2,4)
    pcolor(x,t,psi2)
    shading interp
    colormap jet
    c = colorbar;
    xlabel('x [m]')
    ylabel('t [s]')
    ylabel(c,'|\psi|^2')

    subplot(2,2,3)
    h = plot(x,psi2(1,:),'b.');
    grid
    xlabel('x [m]')
    ylabel('|\psi|^2')
    ht = title('t=0 s');
    ylim([min(psi2(:)),max(psi2(:))])
    for i=2:length(t)
        pause(.01)
        if ~ishandle(h)
            break % Arrete l'animation si la fenetre est fermee
        end
        set(h,'YData',psi2(i,:))
        set(ht,'String',sprintf('t=%0.2f s',t(i)))
    end
end
