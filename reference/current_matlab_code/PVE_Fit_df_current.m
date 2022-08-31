% NANOINDENTATION DATA FIT SCRIPT
% Author : MRI (ECU, Oyen Lab) 07/30/2019
% Edited by DMF (Myers Lab, Columbia) Nov. 2021
% Edited by AMM (Oyen Lab, WashU) July 2022

clear all;
close all;
clc;

%% Indenter Parameters
R=12.7E-3;                                        % Indenter radius(m)
% rampTime=2;                                     % ramp time (s)
holdTime=200;                                    % hold time (s)

m = 1;      %index

%% PIUMA DATA FILE NAME
%Enter location of data (find folder of interest, right click top bar,
%click 'copy address as text'

fnam01= '/Users/anniemascot/Desktop/Oyen Lab SU22/results/indent/20220722_indent';
fnam11='test';

% source_dir = fnam01;
% if ~exist(fullfile(source_dir,'excel'))
%     mkdir (fullfile(source_dir,'excel'))
% end
% dest_dir = strcat(fnam01,'/excel');
% source_files = dir(fullfile(source_dir, '*.csv'));
% for i = 1:length(source_files)
%     data = readmatrix(source_files(i).name);
%     xlswrite(fullfile(dest_dir,source_files(i).name), data);
% end

d = dir;
cd(fnam01);
datadir = fnam01;
k = length(dir);
klen = k-3;
leg = strings(klen,1);
for i = 1:klen
    %     d = dir;
    %      cd(datadir);
    % datadir = fnam01;
    % k = length(dir);
    % klen = k-3;
    % leg = strings(klen,1);
    ii = i+3;
    testname = d(ii).name;
    num2str = string(i);
    leg(i) = num2str + ': ' + testname;
    
    % source_dir = fnam01;
    % dest_dir = '/Users/anniemascot/Desktop/Oyen Lab SU22/results/excel';
    % source_files = testname;
    % for i = 1:length(source_files)
    %     data = csvread(fullfile(source_dir,testname));
    %     xlswrite('20220617_i_15gelatin1.CSV');
    % end
    
    savefig = 'y'; % y or n
    
    %% MATRIX SCAN
    for k=2:1:6       %k represents the x axis
        for l=1:1:6       %l represents the y axis
            
            fnam=sprintf(testname,fnam11,k,l);
            
            %% TEST DATA EXTRACT
            % MATLAB automatically skips empty lines using readtable, so this changes
            % the option to read them
            opt = detectImportOptions(testname);
            opt.EmptyLineRule = 'read';
            % M = readmatrix(testname,'Sheet','Sheet 1','Range','B55:B1000');
            DT = readtable(testname, opt);
            
            % place table variables into a single array called DT
            A = DT.Var2;
            B = DT.Var3;
            C = DT.Var4;
            A = A(52:end);
            B = B(52:end);
            C = C(52:end);
            DT = [A, B, C];
            % Need just first three column data
            DT=[DT(:,1) DT(:,3)*-1 DT(:,2)*1E-3]; % FORMAT- Time(s) Load(N) Indentation(m)
            
            %% DATA ADJUSTMENTS
            % Time Adjust
            %             DT(:,1)=DT(:,1)-DT(1,1);
            %             % Force Adjust
            %             if DT(1,2)<0
            %                 DT(:,2)=DT(:,2)-DT(1,2);
            %             end
            
            %% TEST REGIMES
            % Ftol= 0.001;                                   % force threshold to identify the initial contact point
            % rid=min(find(DT(:,2)>=Ftol))-1;
            index = find(isnan(DT(:,1)));
            contactid = 1;
            dwellid = index(1)+1;
            endid = length(DT);
            % if rid==0 rid=1; end
            [Pmax,pid]=max(abs(DT(:,2)));                  % max load
            rampTime=DT(dwellid,1)-DT(contactid,1);                 % actual ramp time
            ht=DT(endid,1)-DT(dwellid,1); % actual hold time
            %  hid=min(find(DT(:,1)>=ht));
            
            DT1=abs(DT(contactid:1:dwellid-2,:));                      % RAMP  DATA
            DT2=abs(DT(dwellid:1:endid,:));                      % RELAX DATA
            
            % offsetting data to zero
            DT1=DT1-DT1(1,:);
            DT2(:,1)=DT2(:,1)-DT2(1,1);
            DT2(:,3)=DT2(:,3)-DT1(1,3);
            DT2(:,2)=DT2(:,2);                          % unit (N)
            DT2(:,3)=DT2(:,3)*1E-3;                          % unit (m)
            
            % Processed Data Check Plot for Fitting
            % plot(DT2(:,1),DT2(:,2),'-g'); hold on
            
            
            %% PoroViscoelastic Model FIT
            DT2=DT2(1:1:end,:);            % using every data point for fitting
            plotflag=1;
            RR(m,1) = k;    %x position
            RR(m,2) = l;    %y-position
            [RR(m,3:17) PD]=PoroViscoElastic_Model_df_AM(DT2,R,rampTime,holdTime,plotflag,i,leg,dwellid);
            
            m = m+1;
            %Write all results to a file
            T = array2table(RR);
            T.Properties.VariableNames(1:17) = {'X','Y','G0 (kPa)','E0 (kPa)','Ginf (kPa)','Einf (kPa)','Ginf/G0','T1','T2','R-SQ, VISCO','G (kPa)','E (kPa)','nu','k (m^2)','D (m^2/s)','R-SQ, PORO', 'R-SQ, PVE'};
            
            %             if ~exist(fullfile(datadir,'Analysis'))
            %                 mkdir (fullfile(datadir,'Analysis'))
            %             end
            %             save_table = strcat(fnam01,'/Analysis');
            %             cd(save_table)
            if ~exist(fullfile(datadir,'Analysis'))
                mkdir (fullfile(datadir,'Analysis'))
            end
            cd(strcat(fnam01,'/Analysis'))
            writetable(T,strcat(testname, 'results.csv'));
            cd(fnam01)
        end
    end
    
    if savefig == 'y'
        if ~exist(fullfile(datadir,'PVE_fit_plots'))
            mkdir (fullfile(datadir,'PVE_fit_plots'))
        end
        h=findobj('type','figure'); % find the handles of the opened figures
        plotfolder=[datadir '/PVE_fit_plots'];  % Desination folder
        for k=1:numel(h)
            plotfilename=sprintf(strcat(testname,'plot.pdf'),k);
            plotfile=fullfile(plotfolder,plotfilename);
            saveas(h(k),plotfile)
        end
        saveas(h(1),fullfile(plotfolder, 'All.fig'))
    end
end


