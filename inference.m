function figdata=inference(gwas,minv,ns,her2)
  
% Function returns structure figdata with the all the information to plot
% figures. This allows us to create figures with multiple traits if
% needed.

v=gwas.v;
v=v(find(v>minv));
v=sort(v);

disp(['number of included SNPs is',num2str(length(v))]);

disp('regular');
%Getting the main parameters
[vst,U,vtot]=getstat(v,minv);

%For the figures we want the number of variants and their contribution to
%variance as a function of the cutoff contribution to variance.
rng=minv*10.^(-2:0.01:2);
nofv=U*4*expint(2*sqrt(rng/vst));
Vofv=U*vst*exp(-2*sqrt(rng/vst)).*(1+2*sqrt(rng/vst));

nofv0=U*4*expint(2*sqrt(minv/vst));
Vofv0=U*vst*exp(-2*sqrt(minv/vst)).*(1+2*sqrt(minv/vst));

% We calculate the inferred model's CDF to draw parametric bootstrap
% samples from.
vrng=minv*10.^(-2:0.001:4);
cdf0=1-expint(2*sqrt(vrng/vst))./expint(2*sqrt(0.01*minv/vst));
idx=find(diff(cdf0)<=0,1,'first');
if isempty(idx)
    idx=length(cdf0);
end

% We use parametric bootstrap to estimate CI for the inferred model above
% the current v*
boot=100;
for i=1:boot
    if mod(i,1000)==0 
        disp(i);
    end
    v_b=paramboot(vst,minv,U);
    Vofv_b(i,:)=arrayfun(@(mv) sum(v_b.*(v_b>mv)),rng);
end

% We use non-parametric bootsrap to simulate error in parameter estimation
% and parametric bootstrap to estimate sampling noise in order to estimate
% CI for our predictions for future GWAS.
boot=100;
for i=1:boot
    if mod(i,1000)==0 
        disp(i);
    end
    v_b=v(randi(length(v),[length(v),1]));
    [vst_b(i),U_b(i),vtot_b(i)]=getstat(v_b,minv);
    v_b=paramboot(vst_b(i),0.01*minv,U_b(i));

    v_b=v_b(v_b<minv);
    v_b=[v_b;v];
    v_b=sort(v_b);
    nofv_b_over(i,:)=arrayfun(@(mv) sum(v_b>mv),rng);
    Vofv_b_over(i,:)=arrayfun(@(mv) sum(v_b.*(v_b>mv)),rng);
end

% Figure showing the model fit to data
tmp=prctile(Vofv_b,[2.5,97.5]);
idx=find(rng>=minv);
figure;h=plot(rng(idx)/vst,arrayfun(@(mv) 100*sum(v.*(v>mv)),rng(idx))/her2,'r');hold on;h2=plot(rng(idx)/vst,100*Vofv(idx)/her2,'color','k');%[0,0.447,0.741]
xfill=[rng(idx),fliplr(rng(idx))]/vst;yfill=[tmp(1,idx),fliplr(tmp(2,idx))]/her2;
figdata.fig1.xreal=rng(idx)/vst; figdata.fig1.yreal=arrayfun(@(mv) 100*sum(v.*(v>mv)),rng(idx))/her2; figdata.fig1.xsim=rng(idx)/vst;  figdata.fig1.ysim=100*Vofv(idx)/her2; figdata.fig1.xfill=xfill;  figdata.fig1.yfill=yfill; 
hold on;
patch(xfill,100*yfill,1,'FaceColor',[0.8 0.8 0.8],'EdgeColor','none');%[0,0.447,0.741]
%alpha(0.3);
uistack(h2, 'top');
uistack(h, 'top');
set(gca,'layer','top');
xlabel('Threshold variance (relative to v_s)');
%xlabel('Threshold variance');
ylabel('% heritability from loci > v');
%ylabel('% heritability');
standardize(24,3.5,15);
xlim([0 10]);%xlim([0,1e-3*floor(max(v)*1e3)]);

% Figure showing our predictions of the heritability captured in future GWAS
tmp=prctile(Vofv_b_over,[2.5,97.5]);
idxfut=find(rng<=minv);
figure;h=plot(ns*minv./rng(idx),arrayfun(@(mv) 100*sum(v.*(v>mv)),rng(idx))/her2,'r');hold on;h2=plot(ns*minv./rng(idxfut),100*(Vofv(idxfut)-Vofv0+sum(v.*(v>=minv)))/her2,'k');
xfill=[ns*minv./rng(idxfut),fliplr(ns*minv./rng(idxfut))];yfill=[tmp(1,(idxfut)),fliplr(tmp(2,(idxfut)))]/her2;
figdata.fig2.xreal=ns*minv./rng(idx); figdata.fig2.yreal=arrayfun(@(mv) 100*sum(v.*(v>mv)),rng(idx))/her2; figdata.fig2.xsim=ns*minv./rng(idxfut);  figdata.fig2.ysim=100*(Vofv(idxfut)-Vofv0+sum(v.*(v>=minv)))/her2; figdata.fig2.xfill=xfill;  figdata.fig2.yfill=yfill; 
hold on;
patch(xfill,100*yfill,1,'FaceColor',[0.8 0.8 0.8],'EdgeColor','none');
uistack(h2, 'top');
uistack(h, 'top');
set(gca,'layer','top');
xlabel('Study size (in thousands)');
ylabel('Heritability (%)');
standardize(24,3.5,15);
xlim([0,1000]);
set(gca,'xtick',0:200:1000);

% Figure showing our predictions of the number of variants captured in future GWAS
tmp=prctile(nofv_b_over,[2.5,97.5]);
figure;h=plot(ns*minv./rng(idx),arrayfun(@(mv) sum(v>mv),rng(idx)),'r');hold on;h2=plot(ns*minv./rng(idxfut),(nofv(idxfut)-nofv0+sum(v>=minv)),'k');
xfill=[ns*minv./rng(idxfut),fliplr(ns*minv./rng(idxfut))];yfill=[tmp(1,idxfut),fliplr(tmp(2,idxfut))];
figdata.fig3.xreal=ns*minv./rng(idx); figdata.fig3.yreal=arrayfun(@(mv) sum(v>mv),rng(idx)); figdata.fig3.xsim=ns*minv./rng(idxfut);  figdata.fig3.ysim=(nofv(idxfut)-nofv0+sum(v>=minv)); figdata.fig3.xfill=xfill;  figdata.fig3.yfill=yfill; 
hold on;
hp=patch(xfill,yfill+1e-20,1,'FaceColor',[0.8 0.8 0.8],'EdgeColor','none');
uistack(h2, 'top');
uistack(h, 'top');
set(gca,'layer','top');
xlabel('Study size (in thousands)');
ylabel('Number of variants');
set(gca,'yscale','log');
set(gca,'xtick',0:200:1000);
xlim([0,1000]);
standardize(24,3.5,15);

%q-q plots and p_K-S

for dim=[1 10]
    [vst,U,vtot]=getstat(v,minv,dim);
    [quant,quantv]=getquantiles(v,minv,vst,dim);
    
    if 1
    figure; plot(1e4*quantv,1e4*quant);
    hold on; plot([0 100],[0 100],'--k');
    xlim([0 15]);
    ylim([0 15]);
    if dim==1
        title('No pleiotropy');
    else
        title('High pleiotropy');
    end
    standardize(24,2.5,15);
    end
    pks=getpks(v,minv,dim,1000);
    if dim==1
        disp(['The K-S p-value without pleiotropy is',num2str(pks)]);
    else
        disp(['The K-S p-value with pleiotropy is',num2str(pks)]);
    end
    
end
    




% inference with number of traits as free parameter. To use, change 0 to 1. Be warned, it is slow - let run and go
% watch a good movie (better a trilogy).
if 0
[vst,U,vtot,n]=getstatn(v,minv);
disp(vst);
disp(U);
disp(vtot);
disp(n);

n_b=[];
boot=100;
for i=1:boot
    disp(i);
    v_b=v(randi(length(v),[length(v),1]));
    [vst,U,vtot,n]=getstatn(v_b,minv);
    disp(n);
    n_b=[n_b,n];
end
%CI for n
disp(prctile(n_b,[5,100]));
end



end










    


function [pks]=getpks(v,minv,dim,nruns)
disp(minv);
idx=find(v>=minv);
disp(length(idx));
v=sort(v(idx));
nn=length(v);
Dreal=getD(v,minv,dim);
[vst,~,~]=getstat(v,minv,dim);
disp(vst);
vrng=minv*10.^(-1:0.01:3)';
if dim==1
  ecdf=1-expint(2*vrng/vst)./expint(2*0.1*minv/vst);  
else
  ecdf=1-expint(2*sqrt(vrng/vst))./expint(2*sqrt(0.1*minv/vst));
end
idx=find(ecdf<0.99999999);
D=[];vstD=[];vbar=[];
for i=1:nruns
    if mod(i,1000)==0
        disp(i)
    end
    vboot=[];
    while length(vboot)<length(v)
        vtmp=pchip(ecdf(idx),vrng(idx),rand(10*length(v),1));
        idxt=find(vtmp>=minv);
        vboot=[vboot;vtmp(idxt)];
    end
    %disp([length(vboot),length(v)]);
    vboot=vboot(1:length(v));
    %disp([length(vboot),length(v)]);
    %vboot=pchip(ecdf(idx),vrng(idx),rand(size(v)));
    %figure; plot(sort(vboot),sum(vboot)-cumsum(sort(vboot)),'.','markersize',15);
    
    D=[D,getD(vboot,minv,dim)];
    [tmp,~,~]=getstat(vboot,minv,dim);
    vstD=[vstD,tmp];
    if dim==1
      vbar=[vbar,min(vboot)];
    else
      vbar=[vbar,min(sqrt(vboot))];
    end
end
%figure; cdfplot(D); title([num2str(dim),' Dreal=',num2str(Dreal)]);
%figure; cdfplot(vstD); title([num2str(dim),' Dreal=',num2str(vst)]);
%figure; cdfplot(vbar);
pks=length(find(D>Dreal))/length(D);
end

function D=getD(v,minv,dim)
idx=find(v>=minv);
v=sort(v(idx));
nn=length(v);
[vst,~,~]=getstat(v,minv,dim);
if dim==1 
  infcdf=1-expint(2*v/vst)./expint(2*minv/vst);  
else
  infcdf=1-expint(2*sqrt(v/vst))./expint(2*sqrt(minv/vst));
end
D=max(abs([infcdf-(0:nn-1)'/nn;infcdf-(1:nn)'/nn]));%sum(((ecdf-(1:nn)'/nn).^2));%sum(((ecdf-(1:nn)'/nn).^2)./(ecdf.*(1-ecdf)));%max(abs([ecdf-(0:nn-1)'/nn;ecdf-(1:nn)'/nn]));%sum(((ecdf-(1:nn)'/nn).^2)./((1-ecdf)));%max(abs([ecdf-(0:nn-1)'/nn;ecdf-(1:nn)'/nn]));
end

function [quant,quantv]=getquantiles(v,minv,vst,dim)
idx=find(v>=minv);
v=sort(v(idx));
nn=length(v);
vrng=minv*10.^(0:0.01:3)';
if dim==1
  ecdf=1-expint(2*vrng/vst)./expint(2*minv/vst);  
else
  ecdf=1-expint(2*sqrt(vrng/vst))./expint(2*sqrt(minv/vst));
end
idx=find(ecdf<0.9999);
quant=pchip(ecdf(idx),vrng(idx),0.01:0.01:0.99);
quantv=pchip((1:nn)'/nn,v,0.01:0.01:0.99);
end


function [vst,U,vtot]=getstat(v,minv,dim)
    if ~exist('dim','var')
        dim=10;
    end
    if dim==1
        [vst,U,vtot]=getstatConst1D(v,minv);
    else
        [vst,U,vtot]=getstatConst(v,minv);
    end
end



function [vst,U,vtot]=getstatConst(v,minv)
idx=find(v>=minv);
v=sort(v(idx));
sqv=mean(sqrt(v));
%disp(sqv);
opts=optimset('TolX',1e-7);
vst=fminbnd(@(vst) 2*sqv/sqrt(vst)+log(expint(2*sqrt(minv/vst)))     , 0,max(v),opts);
U=length(v)/(4*expint(2*sqrt(minv/vst)));
vtot=U*vst;
end


function [vst,U,vtot]=getstatConst1D(v,minv)
idx=find(v>=minv);
v=sort(v(idx));
meanv=mean(v);
opts=optimset('TolX',1e-7);
vst=fminbnd(@(vst) 2*meanv/vst+log(expint(2*minv/vst))     , 0,max(v),opts);
U=length(v)/(2*expint(2*minv/vst));
vtot=U*vst;
end


function [vst,U,vtot,n]=getstatn(v,minv)
F=@(mvovs,n) 2*integral2(@(v,a) exp(-2*v./(n*a.^2)).*p(a,n)./v,mvovs,Inf,0,1);
f=@(vovs,n) sum(arrayfun(@(vv) log(2*integral(@(a) exp(-2*vv./(n*a.^2)).*p(a,n),0,1)),vovs));


idx=find(v>=minv);
v=sort(v(idx));
sqv=mean(sqrt(v));
%disp(sqv);
opts=optimset('MaxFunEvals',1000,'Display','none');%optimset('MaxFunEvals',1000,'Display','iter');
xin=[1.5;5];
fn=@(x) (abs(x(2))+1)/(0.001*abs(x(2))+1);%  1+ 99*(exp(abs(x))-1)/(1+exp(abs(x)))
fvs=@(x) abs(x(1))*1e-4;
xout=fminsearch(@(x) -f(v/fvs(x),fn(x) )+length(v)*log(F(minv/fvs(x),fn(x) ))     , xin,opts);
vst=fvs(xout);
n=fn(xout);
U=length(v)/(2*F(minv/vst,n));
vtot=U*vst;
end


function v_b=paramboot(vst,minv,U)
persistent vstb
persistent minvb
persistent cdfb
persistent vrng

plotflag=isempty(cdfb);

if isempty(cdfb)||(vstb~=vst)||(minvb~=minv)
    vstb=vst;
    minvb=minv;
    vrng=minvb*10.^(0:0.001:6);
    cdfb=1-expint(2*sqrt(vrng/vstb))./expint(2*sqrt(minvb/vstb));
end

idx=find(diff(cdfb)<=0,1,'first');
if isempty(idx)
    idx=length(cdfb);
end

npois=poissrnd(U*4*expint(2*sqrt(minvb/vstb)));
v_b=pchip(cdfb(1:idx),vrng(1:idx),rand([npois,1]));
end

