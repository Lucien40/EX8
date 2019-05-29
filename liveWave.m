function liveWave(x,t,f,m)
ntsteps = length(t);
dtm = diff(t)/m;
rt = round(t);
figure;
ax = gca;
ax.YLim = [min(min(f)) max(max(f))];
hold on
line = plot(ax,[min(x) max(x)],[0 0],'k:');
line.HandleVisibility = 'off';
for i = 1:ntsteps
    fi = f(i,:);
%     [pks,locs] = findpeaks(fi);
%     [pks2,locs2] = findpeaks(-fi);
    title(['t = ' num2str(rt(i))])
    plot(ax,x,fi,'b.')
%     plot(ax,x(locs),pks,'ro')
%     plot(ax,x(locs2),-pks2,'co')
    if i == ntsteps
        break
    end
    pause(dtm(i)), cla(ax)   
end
hold off
end


