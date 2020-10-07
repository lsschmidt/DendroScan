%%% DendroScan -- Software for multivariate statistical particle shape analyses 
%   after Durig et al. (2012, 2020)
%
%  Copyright (C) 2020 by Tobias Durig and Louise Steffensen Schmidt
%
%  Authors: T. Durig, Institute of Earth Sciences, University of Iceland,
%          Reykjavik, Iceland
%           L. S. Schmidt, Department of Geosciences, University of Oslo,
%           Oslo, Norway
%
%  This program is free software licensed under the GPL (>=v3).
%  Read the GPS-3.TXT file that comes with this software for details.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  This program is free software; you can redistribute it and/or modify it
%  under the terms of the GNU General Public License as published by the
%  Free Software Foundation; either version 3 of the License, or (at your
%  option) any later version.
%  
%  Parts of DendroScan are not copyright by the DendroScan development team.
%  The original authors hold the copyrights and you have to abide to their
%  licensing terms where noted. See the headers of the respective .m scripts
%  for details.
%  
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License (GPL) for more details.
%  
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, see <http://www.gnu.org/licenses/>.
%
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Input is in the form of flat ASCII files in '.csv' format, containing 
%       a list of morphometric parameters. DendroScan was designed to
%       digest tables produced by the software PARTIcle Shape ANalyzer 
%       (PARTISAN)by D�rig et al. (2018).
%      
%  Output files according to selection in the GUI.
%
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  References:
%   Durig et al. (2012): doi.org/10.1007/s00445-011-0562-0
%   Durig et al. (2018): doi.org/10.4401/ag-7865
%   Durig et al. (2020): in review
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If you wish to contribute to the development of DendroScan,
% or to report bugs or other problems with the software, please contact me 
% per email (tobi@hi.is).

function main


warning('off')
mkdir results
warning('on')

addpath functions
hfig=figure(1);
set(hfig,'position',[100 100 1500 500])
movegui(hfig,'center')
start_buttons

end


function start_buttons
clf('reset')
bg = uibuttongroup('Visible','on',...
                  'Position',[0.01 0.6 .1 0.38],...
                  'SelectionChangedFcn',@bselection,'Title','Select test');
              
% Create three radio buttons in the button group.
r(1) = uicontrol(bg,'Style','radiobutton',...
                  'String','none',...
                  'Position',[10 10 100 30],...
                  'HandleVisibility','off');     
              
r(2) = uicontrol(bg,'Style',...
                  'radiobutton',...
                  'String','Dendrogram',...
                  'Position',[10 130 100 30],...
                  'HandleVisibility','off');
              
r(3) = uicontrol(bg,'Style','radiobutton',...
                  'String','T-test',...
                  'Position',[10 100 100 30],...
                  'HandleVisibility','off');

r(4) = uicontrol(bg,'Style','radiobutton',...
                  'String','E-test',...
                  'Position',[10 70 100 30],...
                  'HandleVisibility','on');
r(5) = uicontrol(bg,'Style','radiobutton',...
                  'String','Automatic DAPM',...
                  'Position',[10 40 140 30],...
                  'HandleVisibility','off');

set(figure(1), 'Name', 'DendroScan', 'NumberTitle', 'off')
         
% Make the uibuttongroup visible after creating child objects. 
% bg.Visible = 'on';
              
    function bselection(source,event)
       disp(['Selection: ' event.NewValue.String]);    
       ft=event.NewValue.String;
       if strcmp(ft,'Dendrogram')==1
       selfile_den(ft)
       elseif strcmp(ft,'E-test')==1
       selcal_e(ft)
       elseif strcmp(ft,'T-test')==1
       selfile_t(ft)
       elseif strcmp(ft,'Automatic DAPM')==1
       selfile_aut(ft)
       elseif strcmp(ft,'none')==1
       start_buttons   
       end
    end

end


function selfile_den(ft)
bg = uibuttongroup('Visible','on',...
                  'Position',[0.01 0.02 .1 0.55],...
                  'SelectionChangedFcn',@bselection,'Title','Select files');
    h(1).p = uicontrol(bg,'style','pushbutton','units','pixels',...
                'position',[10,140,120,40],'string','load files',...
                'callback',@p_call);
    h(2).p = uicontrol(bg,'style','pushbutton','units','pixels',...
                'position',[10,70,120,40],'string','load X-matrix',...
                'callback',@p_call2);
    % Pushbutton callback
    function p_call(varargin)
    [file,path] = uigetfile('*.csv','Select three or more dendrogram files','MultiSelect', 'on');

    if(~iscell(file))|| (iscell(file) && size(file,2)<3)
       msg = 'ERROR: Please select at least 3 files for the dendrogram.';
       disp(msg)
       msgbox(msg, 'DendroScan', 'error')
    else
    
    for t=1:length(file)
        f(t).name = file{t};
        f(t).folder = path;
    end
    
    alpha=0.05;
    myGUI2(ft,f,0,0,alpha)
    end
    end

    function p_call2(varargin)
        [file,path] = uigetfile('*.csv','Load X-matrix file','MultiSelect', 'off');
        T=readtable([path '/' file]);
        labels = T.Properties.VariableNames;
        X=table2array(T);
        v_mat = squareform(X);
        tree = linkageold(v_mat,'co');
        subplot(1,20,3:19)
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
        
        endbuttons_den(X,tree,labels)
    end

end

function selcal_e(ft)
    bg = uibuttongroup('Visible','on',...
                      'Position',[0.01 0.25 .1 0.33],...
                      'SelectionChangedFcn',@bselection,'Title','Calibration');
    r(1) = uicontrol(bg,'Style','radiobutton',...
      'String','none',...
      'Position',[10 20 100 30],...
      'HandleVisibility','on');    
    r(2) = uicontrol(bg,'Style','radiobutton',...
      'String','Run calibration',...
      'Position',[10 60 120 30],...
      'HandleVisibility','on');     

    r(3) = uicontrol(bg,'Style',...
      'radiobutton',...
      'String','Import cal file',...
      'Position',[10 100 120 30],...
      'HandleVisibility','on');

    function bselection(source,event)
       disp(['Selection: ' event.NewValue.String]);    
       fj=event.NewValue.String;
       if strcmp(fj,'Import cal file')==1
            selfile_e1(ft,fj)
       else
            selfile_e2(ft,fj)
       end
    end
end

function selfile_e1(ft,fj)

    bg = uibuttongroup('Visible','on',...
                  'Position',[0.01 0.02 .1 0.2],...
                  'SelectionChangedFcn',@bselection,'Title','Select files');
    h.p = uicontrol(bg,'style','pushbutton','units','pixels',...
                'position',[15,20,100,40],'string','load files',...
                'callback',@p_call);
    % Pushbutton callback
    function p_call(varargin)
            [file2,path2] = uigetfile('*.csv','Select calibration file','MultiSelect', 'off');
                fc.name = file2;
                fc.folder = path2;
            [file,path] = uigetfile('*.csv','Select TWO comparison files ','MultiSelect', 'on');
            
            if (~iscell(file)) || (iscell(file) && size(file,2)>2)
               msg = 'ERROR: Please select two comparison files.';
               disp(msg)
               msgbox(msg, 'DendroScan', 'error')
            else
            
            for t=1:length(file)
                f(t).name = file{t};
                f(t).folder = path;
            end
            alpha=0.05;
            myGUI2(ft,f,fc,fj,alpha)
            end
    end
end

function selfile_t(ft)
bg = uibuttongroup('Visible','on',...
                  'Position',[0.01 0.02 .1 0.55],...
                  'SelectionChangedFcn',@bselection,'Title','Select files');
    h.p = uicontrol(bg,'style','pushbutton','units','pixels',...
                'position',[10,100,120,40],'string','load files',...
                'callback',@p_call);
    % Pushbutton callback
    function p_call(varargin)
    [file,path] = uigetfile('*.csv','Select TWO T-test files','MultiSelect', 'on');
    if (~iscell(file)) || (iscell(file) && size(file,2)>2)
       msg = 'ERROR: Please select two t-test files.';
       disp(msg)
       msgbox(msg, 'DendroScan', 'error')
    else
    for t=1:length(file)
        f(t).name = file{t};
        f(t).folder = path;
    end
    alpha=0.05;
    myGUI2(ft,f,0,0,alpha)
    end
    end
end

function selfile_e2(ft,fj)
    bg = uibuttongroup('Visible','on',...
                  'Position',[0.01 0.02 .1 0.2],...
                  'SelectionChangedFcn',@bselection,'Title','Select files');
    h.p = uicontrol(bg,'style','pushbutton','units','pixels',...
                'position',[15,20,100,40],'string','load files',...
                'callback',@p_call);
    % Pushbutton callback
    function p_call(varargin)
            [file2,path2] = uigetfile('*.csv','Select standards','MultiSelect', 'on');
            if(~iscell(file2))
               msg = 'ERROR: Please select more than one standard.';
               disp(msg)
               msgbox(msg, 'DendroScan', 'error')
            else
            
            for t=1:length(file2)
                fc(t).name = file2{t};
                fc(t).folder = path2;
            end
            
            [file,path] = uigetfile('*.csv','Select TWO comparison files ','MultiSelect', 'on');
            if (~iscell(file)) || (iscell(file) && size(file,2)>2)
               msg = 'ERROR: Please select two comparison files.';
               disp(msg)
               msgbox(msg, 'DendroScan', 'error')
            else
            for t=1:length(file)
                f(t).name = file{t};
                f(t).folder = path;
            end
            alpha=0.05;
            myGUI2(ft,f,fc,fj,alpha)
            end
            end
    end

end

function selfile_aut(ft)
bg = uibuttongroup('Visible','on',...
                  'Position',[0.01 0.02 .1 0.55],...
                  'SelectionChangedFcn',@bselection,'Title','Select files');
    h(1).p = uicontrol(bg,'style','pushbutton','units','pixels',...
                'position',[10,100,120,40],'string','load files',...
                'callback',@p_call);
%     h(2).p = uicontrol(bg,'style','pushbutton','units','pixels',...
%                 'position',[10,70,120,40],'string','load X-matrix',...
%                 'callback',@p_call2);
    % Pushbutton callback
    function p_call(varargin)
    [file,path] = uigetfile('*.csv','Select three or more dendrogram files','MultiSelect', 'on');
    if(~iscell(file)) || (iscell(file) && size(file,2)<3)
       msg = 'ERROR: Please select at least 3 files for the dendrogram.';
       disp(msg)
       msgbox(msg, 'DendroScan', 'error')
    else
    for t=1:length(file)
        f(t).name = file{t};
        f(t).folder = path;
    end
    
    txt=[];
    for i=1:length(f)
        txt = [txt f(i).name ' '];
    end
    savename = ['results/DAPM_' datestr(now,'yyyymmdd_HHMM')];
    mkdir(savename)
    yourMsg = 'INFO ';
    fid = fopen([savename '/DAPM_LogFile.txt'], 'a');
    fprintf(fid, '%s\n\n',yourMsg);

    yourMsg = [; 'You have chosen ' num2str(length(f)) ' files: ' txt];
    fid = fopen([savename '/DAPM_LogFile.txt'], 'a');
    if fid == -1
      error('Cannot open log file.');
    end
    fprintf(fid, '%s: %s\n', datestr(now, 0), yourMsg);
    fclose(fid); 

    
    alpha=0.05;
    myGUI2(ft,f,0,0,alpha,savename)
    end
    end

end

function myGUI2(ft,f,fc,fj,alpha,savename)
% Create yes/no checkboxes
bg = uibuttongroup('Visible','on',...
                  'Position',[0.12 0.16 .15 0.82],...
                  'SelectionChangedFcn',@bselection,'Title','Select shape parameters');

names = {'Circ_DL' 'Rec_DL' 'Com_DL' 'Elo_DL' 'Circ_CI' 'AR_CI' ...
   'Con_CI' 'Sol_CI' 'Circ_LL' 'Elo_LL' 'AR_LL' 'Con_LL' ...
   'Sol_LL' 'FF' 'AR_LI' 'Con_LI' 'Sol_LI'  'Circ_SC' 'Rec_SC' ...
   'FF_1' 'AR_F' 'AR_SC' 'Reg'};
j=[1 1 1 1 1 1 1 1 1 1 1 1 0 1 1 0 0 1 0 0 1 1 1];

for i=1:12
h.c(i) = uicontrol(bg,'style','checkbox','units','pixels',...
                'position',[30 355-30*(i-1) 100 30],'string',names{i},'Value',j(i));
end
for i=13:23
h.c(i) = uicontrol(bg,'style','checkbox','units','pixels',...
                'position',[130 355-30*(i-1-12) 100 30],'string',names{i},'Value',j(i));
end
% Create OK pushbutton   
h.p = uicontrol(bg,'style','pushbutton','units','pixels',...
                'position',[80,5,70,20],'string','OK',...
                'callback',@p_call);
    % Pushbutton callback
    function p_call(varargin)
        vals = get(h.c,'Value');
        checked = find([vals{:}]);
        if isempty(checked)
            checked = 'none';
        end 
        
        if strcmp(ft,'Dendrogram')==1
            fc=0; fj=0; 
        alpha_button(f,checked,ft,alpha,fc,fj,names)
        elseif strcmp(ft,'Automatic DAPM')==1
            fc=0; fj=0; 
            txt=[];
            for i=1:length(checked); 
                txt = [txt names{checked(i)} ' ']; 
            end
           yourMsg = [; 'You have chosen ' num2str(length(checked)) ' parameters: ' txt];
           fid = fopen([savename '/DAPM_LogFile.txt'], 'a');
           fprintf(fid, '%s: %s\n', datestr(now, 0), yourMsg);
        alpha_button(f,checked,ft,alpha,fc,fj,names,savename)
        else
        alpha_button(f,checked,ft,alpha,fc,fj,names)
        end
            
        
    end
end

function alpha_button(f,checked,ft,alpha,fc,fj,names,savename)
        bg2 = uibuttongroup('Visible','on',...
                  'Position',[0.12 0.02 .15 0.14],...
                  'SelectionChangedFcn',@bselection,'Title','Change alpha');
        
        hp1 = uicontrol(bg2,'Style', 'pushbutton', 'String', 'OK',...
               'Position', [160 20 30 30], 'Callback', @doit);
        hp2 = uicontrol(bg2,'Style', 'edit', 'String', num2str(alpha),...
               'Position', [50 20 100 30]);
        function doit(~,~)    
            alpha=str2num(get(hp2,'string'));
            
            alpha_button(f,checked,ft,alpha,fc,fj,names,savename)
        end

        if strcmp(ft,'Dendrogram')==1
        [X,tree,labels]=dendrogram_funct(f,alpha,checked,0);
        
        mu = 2.4; s = 4.5; 
        fSPI=@(N)exp(-(N-mu)/s)/(s*(1+exp(-(N-mu)/s))^2); 
        SPI = 100 * fSPI(length(f))/fSPI(2);
        subplot(1,29,7) 
        h=bar(categorical({'SPI'}),[SPI 100-SPI]',1,'stacked','EdgeColor','w');
        if SPI < 1/2*100
            set(h,{'FaceColor'},{'r';'w'});
        elseif SPI < 0.75*100
            set(h,{'FaceColor'},{'y';'w'});
        else
            set(h,{'FaceColor'},{'g';'w'});
        end
        
        endbuttons_den(X,tree,labels)
        elseif strcmp(ft,'Automatic DAPM')==1
        yourMsg = ['You use alpha = ' num2str(alpha)];
        fid = fopen([savename '/DAPM_LogFile.txt'], 'a');
        fprintf(fid, '%s\n\n',yourMsg);
        
        [X,tree,labels,fn,group]=dendrogram_funct(f,alpha,checked,1);
    
       
        Xn = X; Xn(Xn == 0) = NaN; Xn=tril(Xn,-1);
        [ix,iy]=find(isnan(Xn) == 1);
        if ~isempty(ix)
        for i=1:length(ix)
        [~,p]=ttest_main(fn([ix(i) iy(i)]),1,alpha,checked);
        tprob(i)=length(p(p<alpha));
        end
        
        txy=zeros(length(X),1);
        xy =unique([ix; iy]);
        for i=1:length(xy)
        txy(xy(i)) = length(find(xy(i) == [ix;iy]));
        end

        yourMsg = 'WARNINGS ';
        fid = fopen([savename '/DAPM_LogFile.txt'], 'a');
        fprintf(fid, '\n%s\n\n',yourMsg);
        
        if length(find(txy>2)) > 0
            l=find(txy>2);
            for i=1:length(l)
            yourMsg = [; 'WARNING: ' labels{l(i)} ' was used for t-test ' num2str(txy(l(i))) ' times. Type I error may occur'];
            fid = fopen([savename '/DAPM_LogFile.txt'], 'a');            
            fprintf(fid, '%s: %s\n', datestr(now, 0), yourMsg);
            fclose(fid); 
            end
        end
        
        %draw dendrogram
        X(X>0)=X(X>0)+max(tprob);
        for i=1:length(ix)
        X(iy(i),ix(i))=tprob(i);
        X(ix(i),iy(i))=tprob(i);
        end
        
        c=cluster(tree,'cutoff',max(max(X))*0.7,'Criterion','distance');

        tree = linkageold(squareform(X),'co');
        denfig(c,tree,labels,txy)
        if sum(tprob == 0) > 0
        text(max(xlim)*0.05,max(ylim)*1,1,'PRELIMINARY DENDROGRAM AFTER T-TEST','color',[0.5,0.5,0.5],'Fontsize',14)
        end
        else 
            txy = NaN;
        end
        
        fignew = figure('Visible','off');
        set(fignew,'position',[200 200 1000 500])
        H=dendrogram(tree,0,'ColorThreshold','default','Labels',labels);
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
        
        if isnan(txy) | sum(tprob == 0) == 0
            saveas(fignew,[savename '/DAPM_dendrogram.eps'])
            saveas(fignew,[savename '/DAPM_dendrogram.png'])
            saveas(fignew,[savename '/DAPM_dendrogram.tif'])
            close(fignew)
            endbuttons_aut(X,tree,labels);  
        else
            saveas(fignew,[savename '/DAPM_dendrogram_prelim.eps'])
            saveas(fignew,[savename '/DAPM_dendrogram_prelim.png'])
            saveas(fignew,[savename '/DAPM_dendrogram_prelim.tif'])
            close(fignew)
            ebutton_aut(X,tree,labels,alpha,checked,fn,savename,txy)
        end
        elseif strcmp(ft,'T-test')==1
            [pL,pT]=ttest_main(f,1,alpha,checked);
            c=categorical(names(checked));
            c=reordercats(c,names(checked));
            pTn =pT(isfinite(pT));
            in = find(pTn<alpha); 
            
            subplot(1,20,14:18); 
            cla reset
            Tfig(pTn,c,alpha,[0.76 0.0351 0.0981 0.094])
            subplot(1,20,7:11); 
            cla reset
            pLn =pL(isfinite(pL));
            levenefig(pLn,c,alpha,[0.307 0.0351 0.0918 0.094]);
            %text(0.7,-2.8,['Significantly different for ' num2str(length(in)) ' of '  num2str(length(pTn)) ' variables' ])
            a=annotation('textbox', [0.43, 0.03, 0.3, 0.05], 'String', ['Significantly different for ' num2str(length(in)) ' of '  num2str(length(pTn)) ' variables' ]);
            a.BackgroundColor = 'k';
            a.Color = 'y';
            a.HorizontalAlignment = 'c';
            a.FontSize = 11.5;
            a.FontWeight = 'bold';
            
            endbuttons_t(pL,pT,names,checked,alpha)
        elseif strcmp(ft,'E-test')==1
        if strcmp(fj,'Import cal file')==1
           Dmax_cal = readmatrix([fc.folder fc.name ]);
           [Dmax,pL]=etest_main(f,1,alpha,checked);
        else
           Dmax_cal=etest_main(fc,1,alpha,checked);
           [Dmax,pL]=etest_main(f,1,alpha,checked);
        end
        Dmaxn=Dmax(isfinite(Dmax)); Dmax_caln=Dmax_cal(isfinite(Dmax));
        j = find(Dmax_caln < Dmaxn);
        x = 1:length(Dmaxn);
        c=categorical(names(checked));
        c=reordercats(c,names(checked));
        hn=subplot(1,20,14:18); 
        cla(hn);
        etestfig(Dmaxn,Dmax_caln,x,j,names,checked,[0.76 0.0351 0.0981 0.094])

        hn2=subplot(1,20,7:11); 
        cla reset
        pLn =pL(isfinite(pL));
        levenefig(pLn,c,alpha,[0.307 0.0351 0.0918 0.094]);
        
        a=annotation('textbox', [0.43, 0.03, 0.3, 0.05], 'String',['Statistically eqivalent for ' num2str(length(Dmaxn)-length(j)) ' of '  num2str(length(Dmaxn)) ' variables']);
        a.BackgroundColor = 'k';
        a.Color = 'y';
        a.HorizontalAlignment = 'c';
        a.FontSize = 11.5;
        a.FontWeight = 'bold';
%         text(0.7,-2.8,['Statistically eqivalent for ' num2str(length(Dmaxn)-length(j)) ' of '  num2str(length(Dmaxn)) ' variables' ])
        endbuttons_e(Dmax_cal,Dmax,names,checked,pL,alpha)
        end
end


function ebutton_aut(X,tree,labels,alpha,checked,f,savename,txy)
        group = cluster(tree,'cutoff',max(max(X))*0.001,'Criterion','distance');
        ku=unique(group);

% Create yes/no checkboxes
    bg = uibuttongroup('Visible','on',...
                  'Position',[0.88 0.58 .11 0.40],...
                  'SelectionChangedFcn',@bselection,'Title','E-test');

    r(1) = uicontrol(bg,'Style','radiobutton',...
      'String','none',...
      'Position',[10 20 150 30],...
      'HandleVisibility','on');    
    r(2) = uicontrol(bg,'Style','radiobutton',...
      'String','Run calibration (all)',...
      'Position',[10 50 150 30],...
      'HandleVisibility','on');     
    r(3) = uicontrol(bg,'Style',...
      'radiobutton',...
      'String','Import cal file (all)',...
      'Position',[10 80 150 30],...
      'HandleVisibility','on');
  
      r(4) = uicontrol(bg,'Style','radiobutton',...
      'String','Run calibration',...
      'Position',[10 110 150 30],...
      'HandleVisibility','on');     
    r(5) = uicontrol(bg,'Style',...
      'radiobutton',...
      'String','Import cal file',...
      'Position',[10 140 150 30],...
      'HandleVisibility','on');

    function bselection(source,event)
       disp(['Selection: ' event.NewValue.String]);    
       fj=event.NewValue.String;
       yourMsg = 'RESULTS ';
       fid = fopen([savename '/DAPM_LogFile.txt'], 'a');
       fprintf(fid, '\n%s\n',yourMsg);
       if strcmp(fj,'Import cal file')==1
%            tp=NaN(length(ku),1);
           jk=1; 
           nn=[]; 
            for i=1:length(ku)
                if length(find(ku(i)==group)) > 1
                txt=[];
                n = find(ku(i)==group);
                for j=1:length(find(ku(i)==group))
                       txt = [txt ' ' labels{n(j)} ' '];
                end
                [file2,path2] = uigetfile('*.csv',['Select calibration file for ' txt],'MultiSelect', 'off');
                fc.name = file2;
                fc.folder = path2;
                Dmax_cal = readmatrix([fc.folder fc.name ]);
                yourMsg = ['You have chosen the following calibration: ' file2];
                fid = fopen([savename '/DAPM_LogFile.txt'], 'a');
                fprintf(fid, '\n%s: %s\n\n', datestr(now, 0), yourMsg);
                for kk = 1:length(n)
                    for jj = kk:length(n)
                        if kk==jj; continue; end

                        nn=[nn; [n(kk) n(jj)]];
                        [Dmax,pL]=etest_main(f([n(kk) n(jj)]),1,alpha,checked);

                        Dmaxn=Dmax(isfinite(Dmax)); Dmax_caln=Dmax_cal(isfinite(Dmax));

                        tp(jk)=size(find(Dmaxn>Dmax_caln),1); 
                        jk = jk+1;
                    end
                end
                end
            end

            
       elseif strcmp(fj,'Import cal file (all)')==1
            jk=1; 
            [file2,path2] = uigetfile('*.csv',['Select calibration file'],'MultiSelect', 'off');
            fc.name = file2;
            fc.folder = path2;
            yourMsg = ['You have chosen the following calibration: ' file2];
            fid = fopen([savename '/DAPM_LogFile.txt'], 'a');
            fprintf(fid, '\n%s: %s\n\n', datestr(now, 0), yourMsg);
            Dmax_cal = readmatrix([fc.folder fc.name ]);
            nn=[];
            for i=1:length(ku)
                if length(find(ku(i)==group)) > 1
                n = find(ku(i)==group);
                for kk = 1:length(n)
                    for jj = kk:length(n)
                        if kk==jj; continue; end
                        
                        nn=[nn; [n(kk) n(jj)]];
                        [Dmax,pL]=etest_main(f([n(kk) n(jj)]),1,alpha,checked);

                        Dmaxn=Dmax(isfinite(Dmax)); Dmax_caln=Dmax_cal(isfinite(Dmax));

                        tp(jk)=size(find(Dmaxn>Dmax_caln),1); 
                        jk = jk+1;
                    end
                end
                end
            end
       elseif strcmp(fj,'Run calibration')==1
%             tp=NaN(length(ku),1);
            nn=[];
            jk=1; 
            for i=1:length(ku)
                if length(find(ku(i)==group)) > 1
                txt=[];
                n = find(ku(i)==group);
                for j=1:length(find(ku(i)==group))
                       txt = [txt ' ' labels{n(j)} ' '];
                end
                [file2,path2] = uigetfile('*.csv',['Select standards for ' txt],'MultiSelect', 'on');
                if(~iscell(file2))
                   msg = 'ERROR: Please select more than one standard.';
                   disp(msg)
                   msgbox(msg, 'DendroScan', 'error')
                   break
                else
                txt = [];
                for t=1:length(file2)
                    fc(t).name = file2{t}; txt=[txt file2{t} ' '];
                    fc(t).folder = path2;
                end
            
                yourMsg = ['You have chosen ' num2str(length(fc)) ' standards: ' txt ];
                fid = fopen([savename '/DAPM_LogFile.txt'], 'a');
                fprintf(fid, '\n%s: %s\n\n', datestr(now, 0), yourMsg);
                 Dmax_cal=etest_main(fc,1,alpha,checked);
                for kk = 1:length(n)
                    for jj = kk:length(n)
                        if kk==jj; continue; end

                        nn=[nn; [n(kk) n(jj)]];
                        [Dmax,pL]=etest_main(f([n(kk) n(jj)]),1,alpha,checked);

                        Dmaxn=Dmax(isfinite(Dmax)); Dmax_caln=Dmax_cal(isfinite(Dmax));

                        tp(jk)=size(find(Dmaxn>Dmax_caln),1); 
                        jk = jk+1;
                    end
                end
                end
                end
            end
            
       elseif strcmp(fj,'Run calibration (all)')==1
            [file2,path2] = uigetfile('*.csv','Select standards','MultiSelect', 'on');
            jk = 1; txt=[];
            if(~iscell(file2))
               msg = 'ERROR: Please select more than one standard.';
               disp(msg)
               msgbox(msg, 'DendroScan', 'error')
            else
            for t=1:length(file2)
                fc(t).name = file2{t}; txt=[txt file2{t} ' '];
                fc(t).folder = path2;
            end
            
            yourMsg = ['You have chosen ' num2str(length(fc)) ' standards: ' txt ];
            fid = fopen([savename '/DAPM_LogFile.txt'], 'a');
            fprintf(fid, '\n%s: %s\n\n', datestr(now, 0), yourMsg);
            Dmax_cal=etest_main(fc,1,alpha,checked);
            %tp=NaN(length(ku),1);
            nn=[];            
            for i=1:length(ku)
                if length(find(ku(i)==group)) > 1
                n = find(ku(i)==group);
                for kk = 1:length(n)
                    for jj = kk:length(n)
                        if kk==jj; continue; end

                        nn=[nn; [n(kk) n(jj)]];
                        [Dmax,~]=etest_main(f([n(kk) n(jj)]),1,alpha,checked);

                        Dmaxn=Dmax(isfinite(Dmax)); Dmax_caln=Dmax_cal(isfinite(Dmax));

                        tp(jk)=size(find(Dmaxn>Dmax_caln),1); 
                        jk = jk+1;
                    end
                end
                end
            end
            end
       end
       if exist('tp')
       if strcmp(fj,'none') == 0
           X = X + nanmax(tp);
           ns = size(X,1);
           X(1:(ns+1):end) = 0;
           for i=1:size(nn,1)
                   X(nn(i,1),nn(i,2)) = tp(i);
                   X(nn(i,2),nn(i,1)) = tp(i);
                   if tp(i) == 0;
                   yourMsg = [labels{nn(i,1)} ' is statistically equivalent with ' labels{nn(i,2)} ' in ' num2str(length(checked)-tp(i)) ' of ' num2str(length(checked)) ' parameters'];
                   else
                   yourMsg = [labels{nn(i,1)} ' is statistically equivalent with ' labels{nn(i,2)} ' in ' num2str(length(checked)-tp(i)) ' of ' num2str(length(checked)) ' parameters'];
                   end
                   fid = fopen([savename '/DAPM_LogFile.txt'], 'a');
                   fprintf(fid, '%s: %s\n', datestr(now, 0), yourMsg);
           end

           tree = linkageold(squareform(X),'co');
           c=cluster(tree,'cutoff',max(max(X))*0.7,'Criterion','distance');
           denfig(c,tree,labels,txy) 
           
            fignew = figure('Visible','off');
            set(fignew,'position',[200 200 1000 500])
            H=dendrogram(tree,0,'ColorThreshold','default','Labels',labels);
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
            saveas(fignew,[savename '/DAPM_dendrogram_final.eps'])
            saveas(fignew,[savename '/DAPM_dendrogram_final.png'])
            saveas(fignew,[savename '/DAPM_dendrogram_final.tif'])
            close(fignew)
            
            writetable(array2table(X,'VariableNames',labels),[savename '/X_matrix.csv'])
            

            endbuttons_aut(X,tree,labels)
           
       end
       end
       
    end


end


function endbuttons_den(X,tree,labels)
% Create yes/no checkboxes
bg = uibuttongroup('Visible','on',...
                  'Position',[0.9 0.55 .09 0.43],...
                  'SelectionChangedFcn',@bselection,'Title','Export');


% Create OK pushbutton   
tr.a = uicontrol(bg,'style','pushbutton','units','pixels',...
                'position',[15,140,100,40],'string','save X matrix',...
                'callback',@p_call);
tr.b = uicontrol(bg,'style','pushbutton','units','pixels',...
                'position',[15,80,100,40],'string','save tree',...
                'callback',@p_call2);
tr.c = uicontrol(bg,'style','pushbutton','units','pixels',...
    'position',[15,20,100,40],'string','save plot',...
    'callback',@p_call3);
    % Pushbutton callback
    function p_call(varargin)
       [file,path] = uiputfile('Xtest.csv');
       writetable(array2table(X,'VariableNames',labels),[path file])
    end

    function p_call2(varargin)
      [file,path] = uiputfile('dendrogram_tree.csv');
      csvwrite([path file],tree)
      save([path file(1:end-4) '_labels.mat'],'labels')
    end

    function p_call3(varargin)
            fignew = figure('Visible','off');
            set(fignew,'position',[200 200 1000 500])
            D = pdist(X);
%             leafOrder = optimalleaforder(tree,D);
            H=dendrogram(tree,0,'ColorThreshold','default','Labels',labels)%,'Reorder',leafOrder);
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
            [file,path] = uiputfile('dendrogram.png');
            saveas(fignew,[path file])
            close(fignew)
    end

    bg2 = uibuttongroup('Visible','on',...
                      'Position',[0.9 0.02 .09 0.38],...
                      'SelectionChangedFcn',@bselection,'Title','restart');

    tb.a = uicontrol(bg2,'style','pushbutton','units','pixels',...
                    'position',[10,80,100,40],'string','clear all',...
                    'callback',@p_c);


    function p_c(varargin)
        start_buttons
    end

            
end

function endbuttons_e(Dmax_cal,Dmax,names,checked,pL,alpha)
% Create figure
% h.f = figure('units','pixels','position',[200,200,300,50],...
%              'toolbar','none','menu','none');
% Create yes/no checkboxes
bg = uibuttongroup('Visible','on',...
                  'Position',[0.88 0.4 .105 0.58],...
                  'SelectionChangedFcn',@bselection,'Title','Options');


% Create OK pushbutton   
tr.a = uicontrol(bg,'style','pushbutton','units','pixels',...
                'position',[10,230,130,40],'string','save Dmax_cal',...
                'callback',@p_call);
tr.b = uicontrol(bg,'style','pushbutton','units','pixels',...
                'position',[10,160,130,40],'string','save Dmax',...
                'callback',@p_call2);
tr.c = uicontrol(bg,'style','pushbutton','units','pixels',...
    'position',[10,90,130,40],'string','save Levene plot',...
    'callback',@p_call3);
tr.d = uicontrol(bg,'style','pushbutton','units','pixels',...
    'position',[10,20,130,40],'string','save e-test plot',...
    'callback',@p_call4);
    % Pushbutton callback
    function p_call(varargin)
       [file,path] = uiputfile('Dmax_cal.csv');
       csvwrite([path file],Dmax_cal)
    end

    function p_call2(varargin)
      [file,path] = uiputfile('Dmax.csv');
      T=table(names',pL,Dmax,Dmax_cal,'VariableNames',{'parameter', 'Levene test', 'Dmax','Dmax_cal'});
      writetable(T,[path file])
    end

    function p_call4(varargin)
            Dmaxn=Dmax(isfinite(Dmax)); Dmax_caln=Dmax_cal(isfinite(Dmax));
            j = find(Dmax_caln < Dmaxn);
            x = 1:length(Dmaxn);
            fignew = figure('Visible','off');
            etestfig(Dmaxn,Dmax_caln,x,j,names,checked,'best')
            [file,path] = uiputfile('Dmax.png');
            saveas(fignew,[path file])
            close(fignew)
    end

    function p_call3(varargin)
        c=categorical(names(checked));
        c=reordercats(c,names(checked));
        pLn =pL(isfinite(pL));  
        fignew = figure('Visible','off');
        levenefig(pLn,c,alpha,'best')
        [file,path] = uiputfile('levene.png');
        saveas(fignew,[path file])
        close(fignew)
    end

    bg2 = uibuttongroup('Visible','on',...
                      'Position',[0.88 0.02 .105 0.38],...
                      'SelectionChangedFcn',@bselection,'Title','restart');

    tb.a = uicontrol(bg2,'style','pushbutton','units','pixels',...
                    'position',[10,80,130,40],'string','clear all',...
                    'callback',@p_c);


    function p_c(varargin)
        start_buttons
    end

            
end

function endbuttons_t(pL,pT,names,checked,alpha)
% Create yes/no checkboxes
bg = uibuttongroup('Visible','on',...
                  'Position',[0.88 0.5 .105 0.48],...
                  'SelectionChangedFcn',@bselection,'Title','Options');


% Create OK pushbutton   
tr.a = uicontrol(bg,'style','pushbutton','units','pixels',...
                'position',[10,160,130,40],'string','save p-values',...
                'callback',@p_call);
tr.b = uicontrol(bg,'style','pushbutton','units','pixels',...
                'position',[10,90,130,40],'string','save Levene plot',...
                'callback',@p_call2);
tr.c = uicontrol(bg,'style','pushbutton','units','pixels',...
    'position',[10,20,130,40],'string','save t-test plot',...
    'callback',@p_call3);
    % Pushbutton callback
    function p_call(varargin)
       [file,path] = uiputfile('T-test.csv');
       T=table(names',pL,pT,'VariableNames',{'parameter', 'Levene test', 'T-test'});
       writetable(T,[path file])
    end

    function p_call3(varargin)
            fignew = figure('Visible','off');
            c=categorical(names(checked));
            c=reordercats(c,names(checked));
            pTn =pT(isfinite(pT));            
            Tfig(pTn,c,alpha,'best')
            [file,path] = uiputfile('ttest.png');
            saveas(fignew,[path file])
            close(fignew)

    end

    function p_call2(varargin)
        fignew = figure('Visible','off');
        c=categorical(names(checked));
        c=reordercats(c,names(checked));
        pLn =pL(isfinite(pL));            
        levenefig(pLn,c,alpha,'best')
        [file,path] = uiputfile('levene.png');
        saveas(fignew,[path file])
        close(fignew)
    end

    bg2 = uibuttongroup('Visible','on',...
                      'Position',[0.88 0.02 .105 0.38],...
                      'SelectionChangedFcn',@bselection,'Title','Restart');

    tb.a = uicontrol(bg2,'style','pushbutton','units','pixels',...
                    'position',[10,80,130,40],'string','clear all',...
                    'callback',@p_c);


    function p_c(varargin)
        start_buttons
    end

            
end

function endbuttons_aut(X,tree,labels)



% Create OK pushbutton   0.88 0.58 .11 0.40

    bg2 = uibuttongroup('Visible','on',...
                      'Position',[0.88 0.02 .11 0.38],...
                      'SelectionChangedFcn',@bselection,'Title','restart');

    tb.a = uicontrol(bg2,'style','pushbutton','units','pixels',...
                    'position',[20,80,100,40],'string','clear all',...
                    'callback',@p_c);


    function p_c(varargin)
        start_buttons
    end

            
end

function levenefig(pLn,c,alpha,de)
            r=[]; b=[];
            for it = 1:length(pLn)
            hb(it)=barh(c(it),pLn(it)); hold on;
            if pLn(it) < alpha
            set(hb(it),'FaceColor','r');
            r = it;  
            else
            set(hb(it),'FaceColor','b');
            b = it; 
            end
            end
            set(gca,'TickLabelInterpreter','none')
            jk = plot([alpha alpha],[c(1) c(end)],'--k','Linewidth',2);
            title('Levene test')
            if size(b,1) > 0 && size(r,1) > 0
            legend([jk hb(b) hb(r)],'alpha','homogeneous','heterogeneous','Location',de)
            elseif size(b,1) > 0
                legend([jk hb(b)],'alpha','homogeneous','Location',de)
            else
                legend([jk hb(r)],'alpha','heterogeneous','Location',de)
            end
            xlabel('p')
            xlim([0 1])
            scale = 0.1;
            pos = get(gca, 'Position');
            pos(2) = pos(2)+scale*pos(4);
            pos(4) = (1-scale)*pos(4);
            set(gca, 'Position', pos)
end


function Tfig(pTn,c,alpha,de)
            r=[]; b=[];
            for it = 1:length(pTn)
            hb(it)=barh(c(it),pTn(it)); hold on;
            if pTn(it) < alpha
            set(hb(it),'FaceColor','r');
            r = it;  
            else
            set(hb(it),'FaceColor','b');
            b = it; 
            end
            end
            set(gca,'TickLabelInterpreter','none')
            jk = plot([alpha alpha],[c(1) c(end)],'--k','Linewidth',2);
            title('T-test')
            if size(b,1) > 0 && size(r,1) > 0
            legend([jk hb(b) hb(r)],'alpha','difference not verified','significantly different','Location',de)
            elseif size(b,1) > 0
                legend([jk hb(b)],'alpha','difference not verified','Location',de)
            else
                legend([jk hb(r)],'alpha','significantly different','Location',de)
            end
            xlabel('p')
            xlim([0 1])
            scale = 0.1;
            pos = get(gca, 'Position');
            pos(2) = pos(2)+scale*pos(4);
            pos(4) = (1-scale)*pos(4);
            set(gca, 'Position', pos)
            
end

function etestfig(Dmaxn,Dmax_caln,x,j,names,checked,de)

        plot(Dmax_caln,1:length(Dmax_caln),'-k'); hold on
        plot(Dmaxn,x,'*b');
        plot(Dmaxn(j),x(j),'*r');
        set(gca,'TickLabelInterpreter','none')
        set(gca,'ytick',[1:length(Dmaxn)],'yticklabel',names(checked))
        xlim([0 1])

        xlabel('Dmax')
        title('E-test')
        xlim([0 1])

        scale = 0.1;
        pos = get(gca, 'Position');
        pos(2) = pos(2)+scale*pos(4);
        pos(4) = (1-scale)*pos(4);
        set(gca, 'Position', pos)
        
        if ~isempty(j) 
        legend('Dmax calibration','Dmax passed','Dmax failed','Location',de) 
        else
        legend('Dmax calibration','Dmax passed','Location',de)    
        end
        
end

function denfig(c,tree,labels,txy)
        mu = 1; s = 2.5; 
        fSRI=@(N)exp(-(N-mu)/s)/(s*(1+exp(-(N-mu)/s))^2); 
        SRI = 100 * fSRI(max(txy))/fSRI(1);
        subplot(1,29,7) 
        h=bar(categorical({'SRI'}),[SRI 100-SRI]',1,'stacked','EdgeColor','w');
        if SRI < 0.3*100
            set(h,{'FaceColor'},{'r';'w'});
        elseif SRI < 0.75*100
            set(h,{'FaceColor'},{'y';'w'});
        else
            set(h,{'FaceColor'},{'g';'w'});
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

end