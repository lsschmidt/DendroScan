function [pL,pT]=ttest_main(f,het,alpha,w)


for i=1:size(f,2)
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

% Levene test and T-test
[pL,pT]=deal(NaN(length(F),1)); 

for kn=1:length(F)
    if sum(kn == w) == 0; continue; end
    data=F{kn};
    pL(kn)=vartestn_new([data{1}; data{2}],[ones(length(data{1}),1); ones(length(data{2}),1)*2],'TestType','LeveneAbsolute','display','off');
    if het==0 && pL(kn) <= alpha
    else
        if pL(kn) < alpha
            [~,pT(kn)] = ttest2(data{1}, data{2},'Alpha',alpha,'Vartype','unequal');
        else
            [~,pT(kn)] = ttest2(data{1}, data{2},'Alpha',alpha,'Vartype','equal'); 
        end
    end
end

end