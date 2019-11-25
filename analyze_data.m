%% Load data
data=importdata('log_Fspl.txt');
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
