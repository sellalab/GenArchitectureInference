function standardize(fontsize,linewidth,markersize)
if ~exist('fontsize','var')
    fontsize=16;
end
if ~exist('linewidth','var')
    linewidth=0;
end
if ~exist('markersize','var')
    markersize=0;
end

set(gca,'fontsize',fontsize)
set(get(gca,'title'),'fontsize',fontsize);
set(get(gca,'xlabel'),'fontsize',fontsize);
set(get(gca,'ylabel'),'fontsize',fontsize);
set(findobj(gcf,'Type','axes','Tag','legend'),'fontsize',fontsize);
box on;
if linewidth>0
set(get(gca,'children'),'linewidth',linewidth);
end
if markersize>0
set(get(gca,'children'),'markersize',markersize);
end

fontname='arial';%'cambria';%'arial';
set(gca,'fontname',fontname)
set(get(gca,'title'),'fontname',fontname);
set(get(gca,'xlabel'),'fontname',fontname);
set(get(gca,'ylabel'),'fontname',fontname);
set(findobj(gcf,'Type','axes','Tag','legend'),'fontname',fontname);