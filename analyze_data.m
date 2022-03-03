%% Load data
clear;
data=importdata('fort.326');
xi=data(:,1);
xj=data(:,2);
xu=data(:,3);
xv=data(:,4);
xf=data(:,5);
xfu=data(:,6);
xfv=data(:,7);
xfspl=data(:,8);
xfsplu=data(:,9);
xfsplv=data(:,10);
for l=1:length(xi)
    u(xi(l),xj(l))=xu(l);
    v(xi(l),xj(l))=xv(l);
    f(xi(l),xj(l))=xf(l);
    fu(xi(l),xj(l))=xfu(l);
    fv(xi(l),xj(l))=xfv(l);
    fspl(xi(l),xj(l))=xfspl(l);
    fsplu(xi(l),xj(l))=xfsplu(l);
    fsplv(xi(l),xj(l))=xfsplv(l);
end
data=importdata('fort.327');
xi=data(:,1);
xj=data(:,2);
xu=data(:,3);
xv=data(:,4);
xf=data(:,5);
xfu=data(:,6);
xfv=data(:,7);
xfspl=data(:,8);
xfsplu=data(:,9);

%% Make plots
figure('Position',[1 1 1536 512],'Color','white');
subplot(1,2,1);
plot(xu,xf,'k','LineWidth',2); hold on;
plot(xu,xfspl,'r--','LineWidth',2);
xlabel('U'); title('F(u,v)');
legend('F','F\_spl');
subplot(1,2,2);
plot(xu,xfu,'k','LineWidth',2); hold on;
plot(xu,xfsplu,'r--','LineWidth',2);
xlabel('U'); title('dF(u,v)/du');
legend('F','F\_spl');
figure('Position',[1 1 1536 512],'Color','white');
subplot(1,3,1);
pixplot(u,v,f); colormap jet;
xlabel('U'); ylabel('V'); title('F(u,v)');
set(gca,'FontSize',18,'XTick',[-1 0 1],'XTickLabel',{'-1.0','0.0','1.0'},...
    'YTick',[-1 0 1],'YTickLabel',{'-1.0','0.0','1.0'});
subplot(1,3,2);
pixplot(u,v,fspl); colormap jet;
xlabel('U'); ylabel('V'); title('F-Spline');
set(gca,'FontSize',18,'XTick',[-1 0 1],'XTickLabel',{'-1.0','0.0','1.0'},...
    'YTick',[-1 0 1],'YTickLabel',{'-1.0','0.0','1.0'});
subplot(1,3,3);
pixplot(u,v,100.*abs(f-fspl)./abs(f)); colormap jet;
xlabel('U'); ylabel('V'); title('% Difference');
set(gca,'FontSize',18,'XTick',[-1 0 1],'XTickLabel',{'-1.0','0.0','1.0'},...
    'YTick',[-1 0 1],'YTickLabel',{'-1.0','0.0','1.0'});
fmax=max(max(100.*abs(f-fspl)./abs(f)));
text(-0.5,0,['Max. Diff: ' num2str(fmax) '%'],'FontSize',18,'Color','white');

figure('Position',[1 1 1536 512],'Color','white');
subplot(1,3,1);
pixplot(u,v,fu); colormap jet;
xlabel('U'); ylabel('V'); title('dF/du(u,v)');
set(gca,'FontSize',18,'XTick',[-1 0 1],'XTickLabel',{'-1.0','0.0','1.0'},...
    'YTick',[-1 0 1],'YTickLabel',{'-1.0','0.0','1.0'});
subplot(1,3,2);
pixplot(u,v,fsplu); colormap jet;
xlabel('U'); ylabel('V'); title('dF/du-Spline');
set(gca,'FontSize',18,'XTick',[-1 0 1],'XTickLabel',{'-1.0','0.0','1.0'},...
    'YTick',[-1 0 1],'YTickLabel',{'-1.0','0.0','1.0'});
subplot(1,3,3);
pixplot(u,v,100.*abs(fu-fsplu)./abs(fu)); colormap jet;
xlabel('U'); ylabel('V'); title('% Difference');
set(gca,'FontSize',18,'XTick',[-1 0 1],'XTickLabel',{'-1.0','0.0','1.0'},...
    'YTick',[-1 0 1],'YTickLabel',{'-1.0','0.0','1.0'});
fmax=max(max(100.*abs(fu-fsplu)./abs(fu)));
text(-0.5,0,['Max. Diff: ' num2str(fmax) '%'],'FontSize',18,'Color','white');

figure('Position',[1 1 1536 512],'Color','white');
subplot(1,3,1);
pixplot(u,v,fv); colormap jet;
xlabel('U'); ylabel('V'); title('dF/dv(u,v)');
set(gca,'FontSize',18,'XTick',[-1 0 1],'XTickLabel',{'-1.0','0.0','1.0'},...
    'YTick',[-1 0 1],'YTickLabel',{'-1.0','0.0','1.0'});
subplot(1,3,2);
pixplot(u,v,fsplv); colormap jet;
xlabel('U'); ylabel('V'); title('dF/dv-Spline');
set(gca,'FontSize',18,'XTick',[-1 0 1],'XTickLabel',{'-1.0','0.0','1.0'},...
    'YTick',[-1 0 1],'YTickLabel',{'-1.0','0.0','1.0'});
subplot(1,3,3);
pixplot(u,v,100.*abs(fv-fsplv)./abs(fv)); colormap jet;
xlabel('U'); ylabel('V'); title('% Difference');
set(gca,'FontSize',18,'XTick',[-1 0 1],'XTickLabel',{'-1.0','0.0','1.0'},...
    'YTick',[-1 0 1],'YTickLabel',{'-1.0','0.0','1.0'});
fmax=max(max(100.*abs(fv-fsplv)./abs(fv)));
text(-0.5,0,['Max. Diff: ' num2str(fmax) '%'],'FontSize',18,'Color','white');