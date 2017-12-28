function choice_of_minv(gwas,minv)

% minv is v*, the cutoff frequency.
% gwas is a table containing x, the MAF, a, effect sizes, and v, the
% contribution to variance.
% The minv is the minv value we are testing.


v=gwas.v;
v=sort(v);

%vmv is v contitional on v>minv
% We get the value of v_s for minv and the bootstrap CI for it.
vmv=v(v>minv);
vmv=sort(vmv);
[vs0,~,~]=getstat(vmv,minv);
boot=10000;
for i=1:boot
    if mod(i,1000)==0 
        disp(i);
    end
    v_b=vmv(randi(length(vmv),[length(vmv),1]));
    [vs_b(i),~,~]=getstat(v_b,minv);
end


%mv is the value of v* we are trying. For each such value we infer v_s. We
%also count how many sites are above mv.
vss=[];
cnts=[];
for mv=minv*(0.5:0.005:2)
    [vs,~,~]=getstat(v,mv);
    vss=[vss,vs];
    cnts=[cnts,length(find(v>mv))];
end

%plot figure to guide choice of v*. As in Fig. S13. 
figure; h1=plot(minv*(0.5:0.005:2)*1e4,1e4*vss,'-');
hold on;
patch([0 3 3 0],1e4*[prctile(vs_b,2.5) prctile(vs_b,2.5) prctile(vs_b,97.5) prctile(vs_b,97.5)],1,'FaceColor',[0.9 0.9 0.9],'EdgeColor','none');
plot([0 3],1e4*vs0*[1 1],'-','color',[0.5 0.5 0.5]);
xlabel('Cutoff contribution to variance (v*)'); ylabel('Estimated v_s');
xlim(minv*[0.5 2]*1e4);
%ylim([1 3.5]);
disp(['vs=',num2str([vs, prctile(vs_b,[2.5,97.5])])]);
yl=ylim;
h2=plot([minv, minv]*1e4,[0,yl(2)],'--k');
uistack(h1, 'top');uistack(h1, 'top');set(gca,'layer','top');
standardize(24,2.5,15);
end


function [vs,U,vtot]=getstat(v,minv,dim)
    if ~exist('dim','var')
        dim=10;
    end
    if dim==1
        [vs,U,vtot]=getstatConst1D(v,minv);
    else
        [vs,U,vtot]=getstatConst(v,minv);
    end
end



function [vs,U,vtot]=getstatConst(v,minv)
idx=find(v>=minv);
v=sort(v(idx));
sqv=mean(sqrt(v));
%disp(sqv);
opts=optimset('TolX',1e-7);
vs=fminbnd(@(vs) 2*sqv/sqrt(vs)+log(expint(2*sqrt(minv/vs)))     , 0,max(v),opts);
U=length(v)/(4*expint(2*sqrt(minv/vs)));
vtot=U*vs;
end


function [vs,U,vtot]=getstatConst1D(v,minv)
idx=find(v>=minv);
v=sort(v(idx));
meanv=mean(v);
opts=optimset('TolX',1e-7);
vs=fminbnd(@(vs) 2*meanv/vs+log(expint(2*minv/vs))     , 0,max(v),opts);
U=length(v)/(2*expint(2*minv/vs));
vtot=U*vs;
end

