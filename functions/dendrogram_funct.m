function [X,tree,labels,f,group]=dendrogram_funct(f,alpha,w,aut)
% figure;

[X,tree,labels,f] = dendro(f,alpha,w,0);

if aut == 1;

if max(max(X))*0.9 > 0
group=cluster(tree,'cutoff',max(max(X))*0.9,'Criterion','distance');
 
n = 1; 
while n == 1 
[cnew,group,n,X] = groupsort(f,group,alpha,w,0.7,X);
f=f(cnew); X=X(cnew,cnew); labels = labels(cnew);
end
end
tree = linkageold(squareform(X),'co');
end

subplot(1,20,7:19)
H=dendrogram(tree,0,'ColorThreshold','default','Labels',labels);%,'Reorder',leafOrder);
set(gca,'TickLabelInterpreter','none')
set(H,'LineWidth',2);
scale = 0.1;
pos = get(gca, 'Position');
pos(2) = pos(2)+scale*pos(4);
pos(4) = (1-scale)*pos(4);
set(gca, 'Position', pos)
ax=gca; 
xtickangle(ax,90)
ax.YGrid = 'on'; ax.GridLineStyle = ':'; 
ylim([0 ax.YLim(2)+ax.YLim(2)*0.05]);
ylabel('Dissimilarity')



% nested function: levene
function [p,f,fk]=levene(v)

f=[]; fk=[];
t=1; 
for n=1:length(v)
    f=[f; v{n}];
    fk=[fk; ones(length(v{n}),1)*t];  
    t=t+1;
end

p=vartestn_new(f,fk,'TestType','LeveneAbsolute','display','off');
end

function [X,tree,labels,f] = dendro(f,alpha,w,aut)
        
for i=1:size(f,2)
    clearvars T
    try
         T = readtable(fullfile(f(i).folder,f(i).name),'HeaderLines',4,'PreserveVariableNames',1);
    catch
        T = readtable(fullfile(f(i).folder,f(i).name),'HeaderLines',4);
    end

Circ_DL{i}	= table2array(T(:,22));
Rec_DL{i}	= table2array(T(:,23));
Com_DL{i}	= table2array(T(:,24));
Elo_DL{i} = table2array(T(:,25));
Circ_CI{i}	= table2array(T(:,27));
AR_CI{i} = table2array(T(:,28));
Con_CI{i} = table2array(T(:,29));
Sol_CI{i} = table2array(T(:,30));
Circ_LL{i} = table2array(T(:,32));
Elo_LL{i} = table2array(T(:,33));
AR_LL{i} = table2array(T(:,34));
Con_LL{i} = table2array(T(:,36));
Sol_LL{i} = table2array(T(:,37));
FF{i} = table2array(T(:,39));
AR_LI{i} = table2array(T(:,40));
Con_LI{i} = table2array(T(:,41));
Sol_LI{i} = table2array(T(:,42));
Circ_SC{i} = table2array(T(:,44));
Rec_SC{i} = table2array(T(:,45));
FF_1{i} = table2array(T(:,46));
AR_F{i} = table2array(T(:,47));
AR_SC{i} = table2array(T(:,48));
Reg{i} = table2array(T(:,49));
end

F{1}=Circ_DL; F{2}=Rec_DL; F{3}=Com_DL; F{4}=Elo_DL; F{5}=Circ_CI; F{6}=AR_CI;
F{7}=Con_CI; F{8}=Sol_CI; F{9}=Circ_LL; F{10}=Elo_LL; F{11}=AR_LL; F{12}=Con_LL;
F{13}=Sol_LL; F{14}=FF; F{15}=AR_LI; F{16}=Con_LI; F{17}=Sol_LI; F{18}=Circ_SC; 
F{19}=Rec_SC; F{20}=FF_1; F{21}=AR_F; F{22}=AR_SC; F{23}=Reg;

for i =1:length(Circ_CI)
S(i) = mean(Circ_DL{i}(1) + Rec_DL{i}(1));
end

% check for missing numbers
prog = waitbar(0,'Calculating statistics');

for i=1:length(w)
    for j=1:length(F{w(i)})
        
    if sum(isnan(F{w(i)}{j})) > 0
        close(prog)
        a=annotation('textbox', [0.43, 0.5, 0.3, 0.05], 'String', ['Missing values in file ' f(j).name]);
        a.HorizontalAlignment = 'c';
        error(['Missing values in file ' f(j).name]) 
    end
    end 
end

% Levene and 1-way ANOVA with post-hoc
pL=NaN(length(w),1); pA=NaN(length(Circ_DL),length(Circ_DL),length(w));
kk=0; 
for kkn=1:length(F)
    if sum(kkn == w) == 0; continue; end
    prog = waitbar(kkn/(length(F)+2),prog,'Calculating statistics');
    kk = kk+1; 
    [pL(kk),fn,gr]=levene(F{kkn});
    if pL(kk) > alpha
        [~,~, stats] = anova1(fn,gr,'off');
        com = multcompare(stats,'display','off','Ctype','tukey-kramer');
        for i=1:size(com,1)
             pA(com(i,1),com(i,2),kk)=com(i,end);
             pA(com(i,2),com(i,1),kk)=com(i,end);
        end
    else
        [~,p,stats]=games_howell(fn,gr,alpha);
        for i=1:length(p)
            for jk=1:length(p)
                pA(stats.gnames(i),stats.gnames(jk),kk)=p(i,jk);
            end
        end
    end
end


prog = waitbar((kkn+1)/(length(F)+2),prog,'Making X-matrix');

% distance matrix
Y=zeros(length(Circ_DL),length(Circ_DL),length(w));
for i=1:length(Circ_DL)
    for jk=1:length(Circ_DL)
        for kk=1:length(w)
            if pA(i,jk,kk)< 0.05 %&& pA(i,jk,kk) > 0
            Y(i,jk,kk)=log10(1+1/pA(i,jk,kk));
            end
        end
        X(i,jk)=sum(Y(i,jk,:));
    end
end


% dendrogram labels
labels=cell(1,length(Circ_DL));
for i=1:length(Circ_DL)
    labels{i} = f(i).name(1:end-4);
end

if aut == 0
[~,IX]=sort(S); 
X=X(IX,IX); 
labels=labels(IX);
f=f(IX);
end

close(prog)
v_mat = squareform(X);
tree = linkageold(v_mat,'co');
    
end


function [cnew,group,bn,X] = groupsort(f,c,alpha,w,th,X)
cl = unique(c); cnew = []; group = []; gp =1; 
bn=1; 
for i=1:length(cl)
t=find(c==cl(i));


if length(t) < 3
    cnew = [cnew; t];
    group = [group; gp*ones(length(t),1)];
    gp = gp+1; 
else
    xn = 0;
    f1=f(t);
    [Xc,treec,labelsc] = dendro(f(t),alpha,w,1);

    if max(max(Xc))*th == 0 
        cnew = [cnew; t];
        group = [group; gp*ones(length(t),1)];
        gp = gp+1; 
        continue; 
    end
    cc=cluster(treec,'cutoff',max(max(Xc))*th,'Criterion','distance');

    for jk = 1:length(t)
        for ik = 1:length(t)
            X(t(jk),t(ik)) = Xc(jk,ik); 
        end
    end
    
    cnew = [cnew; t];
    group = [group; cc+gp];
    gp = max(group)+1;

end

end
if sum(group - c ) == 0 
    bn = 0;
end
end

end
