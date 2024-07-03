%% 
load ~/Dropbox/DATA/RNA/10X/Mouse-Human_final/01_compile/FINAL/Human.mat
load ~/Dropbox/DATA/RNA/10X/Mouse-Human_final/01_compile/FINAL/Mouse.mat
load ~/Dropbox/MATLAB/FUNCS/viridi.mat

load human_umi.mat
load mouse_umi.mat
load human_genes.mat
load mouse_genes.mat
load human_cells.mat
load mouse_cells.mat

hx = 6; 
mx = 3; 

htime = Human(hx).time; 
mtime = Mouse(mx).time; mtime(mtime > 7) = 0;
mdat = Mouse(mx).Dataset; f = ismember(mdat,'M2'); mtime(f) = 0; 

% remove time stamps for rep 1 because day 10 was not properly timestamped
fh1 = ismember(Human(hx).Dataset,'H1');
htime(fh1) = 0; 

goi_h = {'LIN28A','SOX1','HES1','PAX6','OLIG1','OLIG2','NEUROG2','NEUROG1','LHX3','NKX2-2','SOX9','ISL1','MNX1','SLC18A3'};
goi_m = goi_h;

[hmat,~] = makegoiplot(Human(hx).seurcid,Human(hx).idname,goi_h,hgene,hcell,humi,htime);
[mmat,~] = makegoiplot(Mouse(mx).seurcid,Mouse(mx).idname,goi_m,upper(mgene),mcell,mumi,mtime);

% ho = unique(Human(hx).seurcid);
% mo = unique(Mouse(mx).seurcid); 
ho = [6,8,2,5,3,1,4,0,9,7];
mo = [2,4,6,1,0,3,8,5,7,9];

close all

%% time and heat
figure; 
utim_m = unique(mtime); 
utim_m = utim_m(utim_m > 0);
mtime_mat = NaN(length(mo), length(utim_m));
for i = 1:length(mo)
    m = mo(i);
    f = Mouse(mx).seurcid == m; 
    subplot(length(mo),1,i); 
    temp = mtime(f);
    temp = temp(temp > 0); 
    for j = 1:length(utim_m)
        mtime_mat(i,j) = sum(temp == utim_m(j)); 
    end
    histogram(temp); 
    xlim([3.5,6.5])
    ylabel(sprintf('M%d',mo(i)));
end

figure; 
utim_h = unique(htime); 
utim_h = utim_h(utim_h > 0);
htime_mat = NaN(length(ho), length(utim_h));
for i = 1:length(ho)
    h = ho(i);
    f = Human(hx).seurcid == h; 
    subplot(length(ho),1,i); 
    temp = htime(f);
    temp = temp(temp > 0); 
    for j = 1:length(utim_h)
        htime_mat(i,j) = sum(temp == utim_h(j));
    end
    histogram(temp); 
    xlim([7.5,18.5])
    ylabel(sprintf('H%d',ho(i)))
end

mcl = cell(length(mo),1);
for i = 1:length(mcl)
    mcl{i} = sprintf('M%d',mo(i));
end

hcl = cell(length(ho),1);
for i = 1:length(hcl)
    hcl{i} = sprintf('H%d',ho(i));
end

figure;
heatmap(mcl,goi_m,mmat(:,mo+1),'Colormap',viridi,'GridVisible',true,'ColorLimits',[0.1,0.7],'Title','MOUSE')
figure;
heatmap(hcl,goi_h,hmat(:,ho+1),'Colormap',viridi,'GridVisible',true,'ColorLimits',[0.1,0.6],'Title','HUMAN')

mtime_mat = mtime_mat./sum(mtime_mat,2); 
htime_mat = htime_mat./sum(htime_mat,2);
figure; 
heatmap(utim_m, mcl, mtime_mat, 'Colormap', viridi, 'GridVisible', true, 'Title', 'MOUSE','ColorLimits',[0,0.5])
figure; 
heatmap(utim_h, hcl, htime_mat, 'Colormap', viridi, 'GridVisible', true, 'Title', 'HUMAN','ColorLimits',[0,0.5])

figure; 
mtime_mat_sub = mtime_mat(:,1:3); 
mtime_mat_sub = mtime_mat_sub./sum(mtime_mat_sub,2); 
heatmap(utim_m(1:3), mcl, mtime_mat_sub, 'Colormap', viridi, 'GridVisible', true, 'Title', 'MOUSE rep1','ColorLimits',[0,0.7])


