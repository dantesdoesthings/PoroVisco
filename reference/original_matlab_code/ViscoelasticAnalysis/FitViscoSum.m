function [  ] = FitViscoSum( )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Viscoanalysis code for indentation. 
% Adapted from Matteo Galli's script by Daniel Strange
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%References
% Oyen. Analytical techniques for indentation of viscoelastic materials. Philosophical Magazine (2006) vol. 86 (33-35) pp. 5625-5641
% JM Mattice et al.: Spherical indentation load-relaxation of soft biological tissues, J. Mater. Resl, Vol. 21, No. 8, (2006) 2003-2010
% Oyen. Spherical indentation creep following ramp loading. Journal of Materials Research (2005) vol. 20 (8) pp. 2094-2100
% Oyen and Cook. Load-displacement behavior during sharp indentation of viscous-elastic-plastic materials. Journal of Materials Research (2003) vol. 18 (1) pp. 139-150
% Oyen et al. Indentation responses of time-dependent films on stiff substrates. Journal of Materials Research (2004) vol. 19 (8) pp. 2487-2497
% Oyen and Cook. A practical guide for analysis of nanoindentation data. Journal of the mechanical behaviour of biomedical materials (2008) pp. 12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global options

fprintf('********************** Viscoelastic Analysis routine ************************ \n\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% USER INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INSTRUMENT TYPE
instrument = 'Instron';
%instrument = 'Nanoindenter';
%INDENTER TYPE
%DIMENSIONS HAVE TO BE IN THE OUTPUT UNITS OF THE TEST 
%eg (nm for nanoindentation and mm for instron)
%Sphere
indenter.type='Sphere';
indenter.Radius = 15.9/2;  %mm
%indenter.Radius=238860; %nm

%Berkovich
% indenter.type ='Berkovich';
% indenter.Bangle=70.3;
% indenter.gamma=pi/2.;
% indenter.gamma=1.

%Flat Punch
% indenter.type='Flat Punch';
% indenter.Radius=40.587;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TEST CONTROL TYPE
%control='Load';
control='Displacement';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OUTPUT UNITS
%Units = 'GPa';
%Units = 'MPa';
Units = 'kPa';

%% Number of time constants to fit....

numParam = 2; % number of exponentials in the Creep/Relaxation function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% I/O FILES

%folder where input files are found
datadir = uigetdir; % brings up Folder User Interface when code is run
%datadir = 'C:/MyDocuments/InstronFiles'; % datadir can alternatively be set to the file directory containing the datafiles
%output file with all the results in the current folder of datafiles -
%change name otherwise it'll be over-written
outFile=[datadir '/ViscoElasticAnalysis_2P_fitted.csv'];
comments = 'write a comment here'; % use this to write anything you want into the csv file (can't have commas)


%% Options  %These can generally be ignored!!
% not all of these work for both windows and mac 
options.ReZero = true; %Rezero's data (ie first data point (time = 0, displacement = 0, load = 0)

options.percentRamp = 0; % percent of ramp to fit in load control
options.shortName = false; %shortName doesn't work in windows
options.maxAv = true;  %takes an average of the max load/displacement
options.plotCorrection = false; % this doesn't work very well.
options.fRisePoint = 0.5; % Change this if it doesn't find the rise point correctly. (+ to percent total test time where rise time must be below)
options.fUnloadPoint = 0.5; % same as above but for the unload point
options.optimFunc = 'default'; % this can be used to specify other forms of function to be optomised... eg. could be extended to Poroelasticity
options.optim = optimset();  % optimisation toolbox settings
options.extraInfo = true;    % adds some extra info to the output csv.
options.PlotUnfitReg = true; % plots unfit data
options.PlotRawData = false;  % plots raw data
options.PlotLogX = false;
options.PlotLogY = false;
options.PlotOverview = true; % plots all curves on top of each other

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF USER INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% set suff up 
label.Pmax = 'Pmax';
label.Hmax = 'Hmax';

if options.maxAv
    switch control
        case 'Load'
            label.Pmax = 'PmaxAv';
        case 'Displacement'
            label.Hmax = 'HmaxAv';
    end
end

[m, a, unitFactor] = FitVisco_TestSetup(instrument, indenter, Units);


%% Get file list
if strcmp(instrument, 'Instron')
    ext = 'Spec*.csv';
elseif strcmp(instrument, 'Nanoindenter')
    ext = '*.txt';
end

filelist = fdir(datadir,['/' ext]);

if isempty(filelist)
    disp('** Error... no files correctly selected.. **** Have you selected the correct instrument')
    return
end


%% Process files

for j=length(filelist):-1:1 % doing it in reserve preallocates the variables...
    
    % Get filename (and abreviated filename if selected and put it in 
    % struct with file.name, file.label, file.path....)    
    file(j).name = char(filelist(j));
    file(j).label = FitVisco_fileInfo(file(j).name, datadir);
    fprintf(1, 'file=%s \n',file(j).label);

    % extra data from file and put it in format rawdata
    % Depth	 Load 	Time 
    rawdata{j} = FitVisco_getData(file(j).name, instrument);
    
    % ReZero values (if
    % options.ReZero set..)
    expdata{j} = FitVisco_preProcess(rawdata{j}, instrument, control);
    
    % Find out test params (hmax,pmax, rates, risetime,holdtime etc.)
    [testInfo(j) index] = FitVisco_TestParam(expdata{j}, control, indenter);
    fprintf(1,'%s=%e \n',label.Pmax,testInfo(j).Pmax);
    fprintf(1, '%s=%e \n',label.Hmax, testInfo(j).hmax);
    fprintf(1,'rise time=%e \n',testInfo(j).riseTime);
    fprintf(1,'hold time=%e \n',testInfo(j).holdTime);
    fprintf(1,'rate = %e \n', testInfo(j).rate);
      
      
    % Select fitting data (Ramp can only be included in Load Control)
    % fit data is in format [time depednent variable(load/disp)
    if strcmp(control, 'Load')
        startPoint = find(expdata{j}(:,3)>=(1 - options.percentRamp)*testInfo(j).riseTime, 1, 'first');
        fitdata = expdata{j}(startPoint:index.unload, [3 1]);
    else
        startPoint = index.hold;
        fitdata = expdata{j}(index.hold:index.unload, [3 2]);
    end
      
    % Determine starting parameters for optomisation (guesses and bounds)
    [Xguess, LB, UB] = FitVisco_OptParam(numParam, fitdata, m, a,  testInfo(j), control);
             
    % Perform fitting routine
    [results(j,:) fit] = FitVisco_ParamIdent(numParam, fitdata, control, testInfo(j), m, a, Xguess, LB, UB) ;
    
    % Determine goodness of fit
    RSquare(j) = cofDet(fitdata(:,2), fit');
    
    if strcmp(control,'Load') && options.PlotUnfitReg
        [UnFitRegion] = calcUnFitRegion(results(j,:), index, startPoint, expdata{j}, numParam, testInfo(j), a,m);
    else
        UnFitRegion = zeros(0,2); % create dummy matrix which won't be plotted
    end
    
    % Plot against
    
    FitVisco_Plotter(rawdata{j},expdata{j},fitdata,fit,UnFitRegion,control, instrument, file(j).label, j);
    
end


%% Display and Process Results

% Plot overview of stuff if selected
if options.PlotOverview
    FitVisco_PlotOverview(expdata, control, instrument, file)
end

% Output average + std results to screen.

aveResults= mean(results,1);
stdResults= std(results,0,1);

FitVisco_OutputResults(numParam, Units, aveResults, unitFactor, stdResults);

% Write all results to a file
FitVisco_OutputFile(results, aveResults, stdResults, Units, unitFactor, label, control,  instrument, indenter, comments, numParam, outFile, testInfo, file, RSquare);


end


%% Function to setup units and indenter parameters

function [m, a, unitFactor] = FitVisco_TestSetup(instrument, indenter, Units) 


    %factor which multiplies displacement relaxation curve
    switch indenter.type
        case 'Berkovich'
            m=2.;
            a=2 * pi*tand(indenter.Bangle)/2./indenter.gamma^2; %factor of 2's included to make viscoanalysis easier
        case 'Sphere'
            m=3./2.;
            a= 2* 4.*sqrt(indenter.Radius)/3.;
        case 'Flat Punch'
            m=1.;
            a=2* 2*indenter.Radius;
        otherwise 
            fprintf(2, 'Invalid indenter: %s \n', indenter.type);
            return
    end
    
    %% CALCULATE UNIT FACTOR
    %FACTOR TO CONVERT G and Ginf to output Units
    
    % calculate conversion to MPa
    if strcmp(instrument,'Nanoindenter')
        %1e6 to convert nm/(µN)^2  to MPa
        unitFactor=1e6;
    elseif strcmp(instrument, 'Instron')
        %1 to convert mm/N^2 to convert to MPa
        unitFactor=1;
    end
    
    % calculate MPa to chosen units
    switch Units
        case 'GPa'
            unitFactor = unitFactor/1E3;
        case 'MPa'
            unitFactor = unitFactor; %#ok<ASGSL>
        case 'kPa'
            unitFactor = unitFactor*1E3;
        otherwise
            fprintf(2, 'Units error: %s \n', Units)
    end
      
end

%% Function to determine the file label

function [fileLabel, filePath] = FitVisco_fileInfo(filename,datadir)
    
    global options

    % Determine fileLabel
    if options.shortName
        %doesn't work when non label directory
    fileInfo = regexp(filename, datadir, 'split'); % splits filename into specified path, and local path
    fileLabel = fileInfo{2}; % filelabel is just the local directory path+name
    else fileLabel = filename;
    end
    
    filePath = '';
    
end



%% Function to get raw data from files in Disp, Load, Time format
function [rawdata] = FitVisco_getData(filename, instrument)
       
   %nano rawdata has 3 columns Depth (nm)	Load (µN)	Time (s)
   %instron rawdata has columns Depth (mm)  Load(N)     Time (s)
    %extra data from files
    
    if strcmp(instrument, 'Nanoindenter')
        rawdata=dlmread(filename,'',5,0);% this might need to be 8 if its wrong
    elseif strcmp(instrument, 'Instron')
        csvdata=dlmread(filename,',',2,0); % csvdata comes in format of time, - disp, - load
        rawdata=[-csvdata(:,2) -csvdata(:,3) csvdata(:,1)];   
    else
        fprintf(2, 'Invalid Instrument: %s', instrument);
        return
    end

end

%% Function to preprocess and rezero data

function [expdata] = FitVisco_preProcess(rawdata, instrument, control)

global options
    
if options.ReZero

    %zero time: beginning of the test    
    if  strcmp(instrument,'Nanoindenter')
        zero = 1;        
        zeroDisp = rawdata(zero,1);
        zeroLoad = rawdata(zero,2);
        zeroTime = rawdata(zero,3);
        expdata = rawdata - ones(size(rawdata,1),1) * [zeroDisp, zeroLoad, zeroTime];
        
    elseif strcmp(instrument,'Instron') % this could possibly be changed.. depends on input params...
        zeroDisp = rawdata(1,1);
        zeroLoad = rawdata(1,2); % don't need to rezero load as it should have been zeroed previously
        zeroTime = rawdata(1,3); % don't need to rezero load as it begins when load is zeroed.
        expdata = rawdata - ones(size(rawdata,1),1) * [zeroDisp, zeroLoad, zeroTime];
              
    else
        fprintf(2, ' Displacement not zeroed...')
        return                
    end
    
else
    expdata = rawdata;
    zeroTime = find(expdata(:,3) == 0, 1, 'first');
end

fprintf(1,'zero time=%e \n',zeroTime);
        
end

%% Function to determine test paramerters, like Pmax, Risepoints etc..

function [testInfo index] = FitVisco_TestParam(expdata, control, indenter) 

global options

testInfo.hmax = 0;      %maximum displacement value
testInfo.Pmax = 0;      %maxiums load value
testInfo.rate = 0;      %loading rate
testInfo.riseTime = 0;  %riseTime
testInfo.holdTime = 0;  %holdTime
testInfo.strain = 0;    %strain
index.zero=0;       %index when the ramp starts
index.hold=0;       %index when the hold time starts
index.unload=0;     %index when unload begins or test ends.

%hmax and pmax

    % Pmax
    
    if (options.maxAv && strcmp(control,'Load'))
        sortP = sort(expdata(:,2));
        testInfo.Pmax = mean(sortP((round(0.9*length(sortP)):end)));

    else
        testInfo.Pmax=max(expdata(:,2));
    end
    
           
    % Hmax
       
    if (options.maxAv && strcmp(control,'Displacement'))
        sortH = sort(expdata(:,1));
        testInfo.hmax = mean(sortH((round(0.9*length(sortH)):end)));
    else
        testInfo.hmax=max(expdata(:,1));
    end
    
    

%zero time: beginning of the test    
    index.zero = find(expdata(:,3) == 0, 1,'first');
                
%rise time    
    for i=1:length(expdata)
        if strcmp(control,'Load') 
        if (expdata(i,2)>0.995*testInfo.Pmax && expdata(i+1,2)<=testInfo.Pmax && expdata(i,3)>0)
            index.hold=i;
            testInfo.riseTime=expdata(i,3);
            break
            
        end
        end
        
        if strcmp(control,'Displacement')
        if (expdata(i,1)>0.995*testInfo.hmax && expdata(i+1,1)<=testInfo.hmax && expdata(i,3)>0)
        	index.hold=i;
            testInfo.riseTime=expdata(i,3);
            break
        end
        end
    end
    
    % If rise time point is wrong!
    if (index.hold > options.fRisePoint*size(expdata,1) && strcmp(control,'Load'))
        disp('****** correcting rise point *******')
        indextemp = index.hold;
        % calculate an average Pmax to reduce effect of spikes
        sortP = sort(expdata(:,2));
        PmaxAv = mean(sortP((round(0.9*length(sortP)):end)));
        % Step down the percRise if the above doesn't work (although it
        % should
        for percRise = 0.995:-0.001:0.96
            for i=1:size(expdata,1)
                if (expdata(i,2)>percRise*PmaxAv && expdata(i+1,2)<=PmaxAv && expdata(i,3)>0)
                    index.hold=i;
                    testInfo.riseTime=expdata(i,3)-0;
                    break
                end
            end
            if index.hold < options.fRisePoint*size(expdata,1)
                if options.plotCorrection
                    figure(101);
                    subplot(2,1,1)
                    plot(expdata(:,3), expdata(:,2), 'r', expdata(indextemp,3), expdata(indextemp,2), 'bx', expdata(index.hold,3), expdata(index.hold,2), 'bx')
                    title('Rise Point Correction')
                end
                break
            end
        end     
    end
    
    if (index.hold > options.fRisePoint*size(expdata,1) && strcmp(control,'Displacement'))
        disp('****** correcting rise point *******')
        indextemp = index.hold;
        
        %end of rise will be at Pmax!
        [PmaxR index.hold] = max(expdata(:,2));
        testInfo.riseTime =expdata(index.hold,3)-0;
                
        if options.plotCorrection
            figure(101);
            subplot(2,1,1)
            plot(expdata(:,3), expdata(:,1), 'r', expdata(indextemp,3), expdata(indextemp,1), 'bx', expdata(index.hold,3), expdata(index.hold,1), 'bx')
            title('Rise Point Correction')
        end
                        
    end
  
%hold time
    for i=1:length(expdata)
    if strcmp(control,'Load')
        if expdata(end,2) < 0.5*testInfo.Pmax % Checks for unloading portion...
            if (expdata(i+1,2)<0.995*testInfo.Pmax && expdata(i+5,2)<0.99*testInfo.Pmax && expdata(i,3)>testInfo.riseTime) % set at 0.9438 which is a bit wierd...
                index.unload=i;
                testInfo.holdTime=expdata(i,3)-testInfo.riseTime;
                break
            end
        else
            if (expdata(i,3)>testInfo.riseTime)
                index.unload=i;
                testInfo.holdTime=max(expdata(:,3))-testInfo.riseTime;
            end
        end
    end
    
    
    if strcmp(control,'Displacement')
        if expdata(end,1) < 0.5*testInfo.hmax
            if (expdata(i+1,1)<0.995*testInfo.hmax && expdata(i+5,1)<0.99*testInfo.hmax && expdata(i,3)>testInfo.riseTime)
                index.unload=i;
                testInfo.holdTime=expdata(i,3)-testInfo.riseTime;
                break
            end
            
        else
            if (expdata(i,3)>testInfo.riseTime)
                index.unload=i;
                testInfo.holdTime=max(expdata(:,3))-testInfo.riseTime;
            end
        end
    end
    end
    
    % If hold time point is wrong!
    if (index.unload < (options.fUnloadPoint)*size(expdata,1) && strcmp(control,'Load'))
        disp('****** correcting hold point *******')
        %disp(sprintf('last %d, 50p point %d', last, 0.5*size(expdata,1)))
        indextemp = index.unload;
        sortP = sort(expdata(:,2));
        PmaxAv = mean(sortP((round(0.9*length(sortP)):end)));
        for percHold = 0.99:-0.001:0.93
            for i=1:size((expdata),1)
                if (expdata(i+1,2)<0.995*PmaxAv && expdata(i+5,2)<percHold*PmaxAv && expdata(i,3)>testInfo.riseTime)
                    index.unload=i;
                    testInfo.holdTime=expdata(i,3)-testInfo.riseTime;
                    break
                end
            end
            if index.unload > (options.fUnloadPoint)*size(expdata,1)
                fprintf('PercHold %d \n', percHold);
                if options.plotCorrection
                    figure(101);
                    subplot(2,1,2)
                    plot(expdata(:,3), expdata(:,2), 'r', expdata(indextemp,3), expdata(indextemp,2), 'bx', expdata(index.unload,3), expdata(index.unload,2), 'bx')
                    title('Hold Point Correction')
                end
                break
            end
        end
    end
    
    % Determine rate
    if strcmp(control,'Load')
        testInfo.rate=testInfo.Pmax/testInfo.riseTime; 
    elseif strcmp(control,'Displacement')
        testInfo.rate=testInfo.hmax/testInfo.riseTime; 
    end
         
    % Determine strain
    switch indenter.type
        case 'Sphere'
            testInfo.strain = 0.2*sqrt(testInfo.hmax/indenter.Radius);
        otherwise
            testInfo.strain = 0; % this isn't really supported yet
    end
end


%% Choose guess's and lower bounds. 
function [Xguess, LB, UB] = FitVisco_OptParam(numParam, fitdata, m, a, testInfo, control) 

    global options

    Pmax = testInfo.Pmax;
    hmax = testInfo.hmax;
    riseTime = testInfo.riseTime;
    holdTime = testInfo.holdTime;
    
    if strcmp(control,'Load')
        C0guess=a*(fitdata(end,2))^(m)/Pmax;
    elseif strcmp(control,'Displacement')
        C0guess=fitdata(end,2)/(hmax^m)/a;
    end
    
    % Choose guesses based on C0guess and rise and holdTimes
    if numParam == 1
        Tguess = 0.5 * holdTime;
        Cguess = 0.1 * C0guess;
    elseif numParam == 2
        Tguess = [riseTime holdTime];
        Cguess = [0.1 * C0guess, 0.1 *C0guess];
    elseif numParam == 3
        Tguess = [riseTime*0.1 riseTime holdTime];
        %Cguess = [0.3 * C0guess, 0.2 *C0guess, 0.1 *C0guess];
        Cguess = [0.1 * C0guess, 0.1 *C0guess, 0.1 *C0guess];
    elseif numParam == 4
        Tguess = [riseTime*0.1 riseTime holdTime/2 holdTime];
        Cguess = [0.1* C0guess, 0.1* C0guess, 0.1 *C0guess, 0.1*C0guess];
    else
        disp('Automatically generating guess... probably not a very good one')
        Tguess = holdTime .* logspace(-numParam/2,0,numParam);
        Cguess = C0guess  .* logspace(-numParam/2,-0.3,numParam);
    end
    
    % create guesses based on the optomisation function chosen
    switch options.optimFunc
        case 'default'
            Xguess = [C0guess Cguess Tguess];
        case 'Csum'
            Csumguess = C0guess - sum(Cguess);
            Xguess=[Csumguess Cguess Tguess];
        case 'Gdiff'
            disp('****** not done yet ******')
        otherwise
            disp('****** unknown opt function ******')
    end
        
    % Use previous results as guess values if keepPrev is selected.
    % Currently only works for selecting the same numparameters or one
    % more.. but the later is not well implented
%     if keepPrev
%         if size(prevresults,2) == 2*numParam +4
%             Xguess=prevresults(j,1:2*numParam +1);
%         elseif size(prevresults,2) ==2*(numParam-1)+4
%             Xguess=[prevresults(j,1:numParm+1), prevresults(j,1)*0.05,  prevresults(j,numParam+2:2*numParam +1), riseTime*0.1];
%         end
%     end
    
    % Set up the bounds on the data % currently done only for 5 but should
    % be set at the top!
    
    % initial guess * is multiplied by Xguess    
    if numParam == 1
        LB = [0 0 0 ];%zeros(1, 2* numParam + 1); % lowerbound
        UB = [];
    elseif numParam == 2
        LB = [0 0 0 0 0];
        UB = [];
    elseif numParam == 3
        LB = [0 0 0 0 0 0 0];
        UB = [inf inf inf inf 500/Tguess(1) 500/Tguess(2) 500/Tguess(3)];%[inf inf inf];
        if strcmp(options.optimFunc, 'Csum')
            LB = [0 0 0 0 0 0 0];
            %LB = [0.25/Csumguess 0 0 0 0 0 0]; % equivalent to saying G0 can't be greater than 2.5 MPa...
        end%[inf, 1, inf, inf, inf]; %[inf 1.5 1.5 4 4];
    elseif numParam == 4
        LB = [0 0 0 0 0 0 0 0 0];
        UB = [inf inf inf inf inf inf/Tguess(1) inf/Tguess(2) inf/Tguess(3) inf/Tguess(4)];
    else
        LB = zeros(1,2*numParam+1); 
        UB = [];
    end
    
end

%% Perform fitting routine, using fit data and other parameters!


function [results fit] = FitVisco_ParamIdent(numParam, fitdata, control, testInfo, m, a, Xguess, LB, UB)

global options

X0 = ones(1, 2 * numParam + 1); % initial arguments for all functions start at one.

    if strcmp(control,'Load')
                         

        weights = ones(length(fitdata(:,1)),1);

        %fit algorithm
        switch options.optimFunc
        case 'default'
            X = lsqnonlin(@OBJVisLoadC3P,X0,LB,UB,options.optim,fitdata,Xguess,testInfo.riseTime,testInfo.rate,a,m, weights);
            Xid=real(X.*Xguess);
            C0 = Xid(1);
            C=Xid(2:numParam +1);
            T=Xid(numParam + 2:2*numParam + 1);
        case 'Csum'
            X = lsqnonlin(@OBJVisLoadC3P_Csum,X0,LB,UB,options.optim,fitdata,Xguess,testInfo.riseTime,testInfo.rate,a,m, weights);
            Xid=real(X.*Xguess);
            Csum = Xid(1);
            C=Xid(2:numParam +1);
            C0 = Csum+ sum(C);
            T=Xid(numParam + 2:2*numParam + 1);
        case 'Gdiff'
            disp('****** not done yet ******')
        otherwise
            disp('****** unknown opt function ******')
            return
        end
        
        % calculate fit based on identified parameters and plot.
                
        for k=length(fitdata(:,1)):-1:1 % prealocate variable sizes.
        if fitdata(k,1)<testInfo.riseTime
            fit(k)=(1/a)^(1/m)*(C0*testInfo.rate*fitdata(k,1)-sum(C.*T.*testInfo.rate.*(1-exp(-fitdata(k,1)./T))))^(1/m);
        end
        if fitdata(k,1)>=testInfo.riseTime
            fit(k)=(1/a)^(1/m)*(C0*testInfo.rate*testInfo.riseTime-sum(C.*T.*testInfo.rate.*exp(-fitdata(k,1)./T).*(exp(testInfo.riseTime./T)-1)))^(1/m);
        end
        end
        
        % calculate parameters, store and plot
        G0=1/2./(C0-sum(C));
        Ginf=1/2./C0;

        
       
    elseif strcmp(control,'Displacement')
        % fit algorithm - only fits hold portion
         %fit algorithm
        switch options.optimFunc
        case 'default'
            X = lsqnonlin(@OBJVisDispC3P,X0,LB,UB,options.optim,fitdata,Xguess,testInfo.riseTime,testInfo.hmax,a,m);
            Xid = real(X.*Xguess);            
            C0 = Xid(1);
            C=Xid(2:numParam +1);
            T=Xid(numParam + 2:2*numParam + 1);
        case 'Csum'
            X = lsqnonlin(@OBJVisDispC3P_Csum,X0,LB,UB,options.optim,fitdata,Xguess,testInfo.riseTime,testInfo.hmax,a,m);
            Xid = real(X.*Xguess);    
            Csum = Xid(1);
            C=Xid(2:numParam +1);
            C0 = Csum+ sum(C);
            T=Xid(numParam + 2:2*numParam + 1);
        case 'Gdiff'
            disp('****** not done yet ******')
        otherwise
            disp('****** unknown opt function ******')
            return
        end
                
        for k=length(fitdata(:,1)):-1:1
        fit(k)=C0*a*testInfo.hmax^m+sum(C.*a*testInfo.hmax^m.*T./testInfo.riseTime.*(exp(testInfo.riseTime./T)-1).*exp(-fitdata(k,1)./T));
        end
        
        G0=(C0+sum(C))/2.;
        Ginf=C0/2.;
    
    else fprintf(2, 'Unknown control %s \n', control')
        
    end
    
    results= [C0 C T G0 Ginf 4*G0 Ginf/G0]; 
    
    
end

%% Determine unfit portion for load control

function [UnFitRegion] = calcUnFitRegion(results, index, startPoint, expdata, numParam, testInfo, a,m)
 % calculate fit for region that was not fitted based on identified parameters eg ramp?. % this should
 % really be an embbed function....
 
 %Get fit variables
 C0 = results(1);
 C =  results(2:numParam+1);
 T =  results(numParam+2:2*numParam+1);
 
 % select ramp portion time that hasn't yet been fitted...
 UnFitRegion(:,1) = expdata(index.zero:startPoint, 3);
 
  for k=length(UnFitRegion(:,1)):-1:1 % prealocate variable sizes.
     if  k < index.hold - index.zero
         UnFitRegion(k,2)=(1/a)^(1/m)*(C0*testInfo.rate*UnFitRegion(k,1)-sum(C.*T.*testInfo.rate.*(1-exp(-UnFitRegion(k,1)./T))))^(1/m);
     end
     if k >= index.hold - index.zero
         UnFitRegion(k,2)=(1/a)^(1/m)*(C0*testInfo.rate*testInfo.riseTime-sum(C.*T.*testInfo.rate.*exp(-UnFitRegion(k,1)./T).*(exp(testInfo.riseTime./T)-1)))^(1/m);
     end
 end

end


%% Plots results in nice plot...

function [] = FitVisco_Plotter(rawdata,expdata,fitdata, fit, UnFitRegion, control, instrument, filelabel, fignum)
 
global options

%% plot results and fit...

% style parameters
rawDataStyle = 'b'; 
expDataStyle  = 'ko';
fitDataStyle  = 'r';
unFitRegStyle = '--r';
markerSize = 2;
lineThickness = 1.5;

if ~options.PlotRawData
    rawdata = zeros(0,3); % set rawdata to empty matrix so its not plotted
end

% determine units for axis labels.
switch instrument
    case 'Instron'
        load_units = 'N';
        disp_units = 'mm';
    case 'Nanoindenter'
        load_units = '$$\mu$$N';
        disp_units = 'nm';
end

% plot the graph

    figure(fignum); clf;

    if strcmp(control,'Displacement')
        yAxisLabel = sprintf('Load, %s (%s)', '\emph{P}', load_units);
        plotHandles = plot(rawdata(:,3),rawdata(:,2), rawDataStyle, expdata(:,3),expdata(:,2), expDataStyle, fitdata(:,1),fit, fitDataStyle,'MarkerSize',markerSize,'LineWidth',lineThickness);
    elseif strcmp(control, 'Load')
        yAxisLabel = sprintf('Displacement, %s (%s)', '\emph{h}', disp_units);
        plotHandles = plot(rawdata(:,3),rawdata(:,1), rawDataStyle, expdata(:,3),expdata(:,1), expDataStyle, fitdata(:,1),fit, fitDataStyle, UnFitRegion(:,1),UnFitRegion(:,2), unFitRegStyle,'MarkerSize',markerSize,'LineWidth',lineThickness);
    else
        fprintf(2, 'Invalid control %s \n', control);
    end
plotTitle  = title(filelabel,   'fontsize',14, 'FontName'   , 'Helvetica', 'interpreter','none');
plotXlabel = xlabel('Time, \emph{t} (s)',   'fontsize',14,  'interpreter','latex');
plotYLabel = ylabel(yAxisLabel, 'fontsize',14,  'interpreter','latex');

set( gca ,'FontName'   , 'Helvetica' ,'fontsize',11);
set(plotHandles(1) , 'MarkerSize', 2,'LineWidth', 1.25)


if options.PlotLogX
    set(gca,'Xscale','log');
end
if options.PlotLogY
    set(gca,'Yscale','log');
end

end
%% plot exdata

function [] = FitVisco_PlotOverview(expdata, control, instrument, file)

% determine units for axis labels.
switch instrument
    case 'Instron'
        load_units = 'N';
        disp_units = 'mm';
    case 'Nanoindenter'
        load_units = '$$\mu$$N';
        disp_units = 'nm';
end

% clear figure
figure(1000); clf

% begin plotting...
for j = 1:length(expdata)
    
    hold all
    if strcmp(control,'Displacement')
        yAxisLabel = sprintf('Load, %s (%s)', '\emph{P}', load_units);
        plot(expdata{j}(:,3),expdata{j}(:,2),'DisplayName',file(j).label)
    elseif strcmp(control, 'Load')
        yAxisLabel = sprintf('Displacement, %s (%s)', '\emph{h}', disp_units);
        plot(expdata{j}(:,3),expdata{j}(:,1),'DisplayName',file(j).label)
    else
        fprintf(2, 'Invalid control %s \n', control);
    end
    
   hold off
    
end

plotTitle  = title('Plot Overview',   'fontsize',14, 'interpreter','latex');
plotXlabel = xlabel('Time, \emph{t} (s)',   'fontsize',14,  'interpreter','latex');
plotYLabel = ylabel(yAxisLabel, 'fontsize',14,  'interpreter','latex');

plot2leg = legend('Show','Location','SouthEast');
set(plot2leg, 'interpreter','none');

end

%% OUTPUT RESULTS to display screen

function [] = FitVisco_OutputResults(numParam, Units, aveResults, unitFactor, stdResults)


if numParam == 1
    fprintf('Instaneous Shear Modulus G0[%s]: %e +- %e \n',Units, aveResults(4)*unitFactor,stdResults(4)*unitFactor);
    fprintf('Long-term Shear Modulus Ginf[%s]: %e +- %e \n',Units, aveResults(5)*unitFactor,stdResults(5)*unitFactor);
    fprintf('Time constant tau1[s]: %e +- %e \n',aveResults(3),stdResults(3));
    
elseif numParam == 2
    fprintf('Instaneous Shear Modulus G0[%s]: %e +- %e \n',Units, aveResults(6)*unitFactor,stdResults(6)*unitFactor);
    fprintf('Long-term Shear Modulus Ginf[%s]: %e +- %e \n',Units, aveResults(7)*unitFactor,stdResults(7)*unitFactor);
    fprintf('Time constant tau1[s]: %e +- %e \n',aveResults(4),stdResults(4));
    fprintf('Time constant tau2[s]: %e +- %e \n',aveResults(5),stdResults(5));
    
else
    fprintf('Instaneous Shear Modulus G0[%s]: %e +- %e \n',Units, aveResults(end-3)*unitFactor,stdResults(end-3)*unitFactor);
    fprintf('Long-term Shear Modulus Ginf[%s]: %e +- %e \n',Units, aveResults(end-2)*unitFactor,stdResults(end-2)*unitFactor);
    fprintf('Viscoelastic Ratio Ginf/G0: %f +- %f \n', aveResults(end),stdResults(end));
    
    
end

end

%% OUTPUT RESULTS to file\

function [] = FitVisco_OutputFile(results, aveResults, stdResults, Units, unitFactor, label, control, instrument, indenter, comments, numParam, outFile, testInfo, file, RSquare)

global options

%output file
% csv file
fid=fopen(outFile, 'wt');

% CREATE STRINGS TO WRITE TO FILE

% Add extra info to the csv file if requested in options (figurer numbers,
% max depths.. etc.
if options.extraInfo
    extra.Header     = sprintf('%s, %s, %s, %s, %s, %s, ', 'FigNum', 'RSquare', label.Pmax        ,  label.Hmax,          'Strain'           , 'RiseTime');
    extra.AveLine    = sprintf('%s, %.4f, %e, %e, %f, %.2f,', '-'     , mean(RSquare), mean([testInfo.Pmax]), mean([testInfo.hmax]),   mean([testInfo.strain]), mean([testInfo.riseTime])); 
    extra.StdLine    = sprintf('%s, %.4f, %e, %e, %f, %.2f,', '-'     , std(RSquare), std([testInfo.Pmax]) ,  std([testInfo.hmax]) ,   std([testInfo.strain]) , std([testInfo.riseTime]));  
    for i=1:size(results,1)
    extra.results{i} = sprintf('%f, %.4f, %e, %e, %f, %.2f,', i       , RSquare(i), testInfo(i).Pmax   ,  testInfo(i).hmax   ,   testInfo(i).strain,    testInfo(i).riseTime  ); 
    end
else
    extra.Header = '';
    extra.AveLine = '';
    extra.StdLine = '';
    for i=1:size(results,1)
    extra.results{i} = '';
    end
end


% Add model parameters = C0 etc...
if numParam >=1
    if strcmp(control, 'Load')
        CFactor = 1/unitFactor;
        CUnitLabel = sprintf('1/%s', Units);
    elseif strcmp(control,'Displacement')
        CFactor = unitFactor;
        CUnitLabel = Units;
    end
    
    paramWrite.Header  = sprintf('C0[%s], ', CUnitLabel);
    paramWrite.AveLine = sprintf('%e, ', aveResults(1)*CFactor);
    paramWrite.StdLine = sprintf('%e, ', stdResults(1)*CFactor);
    for i=1:size(results,1)
        paramWrite.results{i} = sprintf('%e, ', results(i,1)*CFactor);
    end
    for j = 1:numParam;
        tempWrite.Header  =    sprintf('C%d[%s], T%d[s], ', j, CUnitLabel, j);
        tempWrite.AveLine =    sprintf('%e, %f, ', aveResults(j+1)*CFactor, aveResults(j+1+numParam));
        tempWrite.StdLine =    sprintf('%e, %f, ', stdResults(j+1)*CFactor, stdResults(j+1+numParam));
        for i=1:size(results,1)
            tempWrite.results{i} = sprintf('%e, %f, ', results(i,j+1)*CFactor,   results(i,j+1+numParam));
        end
        
        paramWrite.Header = [paramWrite.Header tempWrite.Header];
        paramWrite.AveLine = [paramWrite.AveLine tempWrite.AveLine];
        paramWrite.StdLine = [paramWrite.StdLine tempWrite.StdLine];
        for i=1:size(results,1)
            paramWrite.results{i} = [paramWrite.results{i} tempWrite.results{i}];
        end
    end
else
    paramWrite.Header = '';
    paramWrite.AveLine = '';
    paramWrite.StdLine = '';
    for i=1:size(results,1)
    paramWrite.results{i} = '';
    end
end


% Create w/ Shear Moduli, reduced moduli, etc...


shearWrite.Header     = sprintf('%s, %s, %s, %s, ',['G0[' Units ']' ]          ,['Ginf[' Units ']' ]        ,['Er[' Units ']' ]          ,'Ginf/G0'    );
shearWrite.AveLine    = sprintf('%f, %f, %f, %f, ',aveResults(end-3)*unitFactor,aveResults(end-2)*unitFactor,aveResults(end-1)*unitFactor,aveResults(end));
shearWrite.StdLine    = sprintf('%f, %f, %f, %f, ',stdResults(end-3)*unitFactor,stdResults(end-2)*unitFactor,stdResults(end-1)*unitFactor,stdResults(end));
for i=1:size(results,1)
shearWrite.results{i} = sprintf('%f, %f, %f, %f, ',results(i,end-3)*unitFactor ,results(i,end-2)*unitFactor ,results(i,end-1)*unitFactor ,results(i,end) );
end

% Options info..... 
optionsInfo{1} = sprintf('Options');
optionsInfo{2} = sprintf('Instrument, %s', instrument);
optionsInfo{3} = sprintf('Tip Type, %s'  , indenter.type);
optionsInfo{4} = sprintf('Indenter Radius, %f ', indenter.Radius); %this won't work for a berk
optionsInfo{5} = sprintf('Control, %s ', control');
optionsInfo{6} = sprintf('Comments:, %s', comments');

% WRITE STRINGS INTO CSV FILE

fprintf(fid,'%s %s %s %s\n' ,extra.Header  ,   paramWrite.Header,     shearWrite.Header,    'filename');
fprintf(fid,'%s %s %s %s\n' ,extra.AveLine ,   paramWrite.AveLine,    shearWrite.AveLine,   'average');
fprintf(fid,'%s %s %s %s\n' ,char(extra.StdLine) ,   char(paramWrite.StdLine),    char(shearWrite.StdLine),   'std dev');
for i=1:size(results,1)
fprintf(fid, '%s %s %s %s\n',extra.results{i}, paramWrite.results{i}, shearWrite.results{i}, file(i).label);  %char(procfiles(i)));
end
% Print some blank lines
for i = 1:5
fprintf(fid, '\n');
end
% add in options info
for i = 1:size(optionsInfo,2)
fprintf(fid, ', , , %s \n', optionsInfo{i});
end

fclose(fid);

% text file
% fid=fopen(outFile, 'wt');
% fprintf(fid,'%s\t\t%s\t\t%s\t\t\t%s\t\t\t%s\t\t%s\t\t%s\t\t%s\n','C0[1/MPa]','C1[1/MPa]','T1[s]','G0[MPa]','Ginf[MPa]','Er[MPa]','Ginf/G0','filename');
% fprintf(fid,'%e\t%e\t%e\t%e\t%e\t%e\t%e\t%s\n',aveResults(1)/unitFactor,aveResults(2)/unitFactor,aveResults(3),aveResults(4)*unitFactor,aveResults(5)*unitFactor,aveResults(6)*unitFactor,aveResults(7),'average');
% fprintf(fid,'%e\t%e\t%e\t%e\t%e\t%e\t%e\t%s\n',stdResults(1)/unitFactor,stdResults(2)/unitFactor,stdResults(3),stdResults(4)*unitFactor,stdResults(5)*unitFactor,stdResults(6)*unitFactor,stdResults(7),'std dev');
% for i=1:size(results,1)
% fprintf(fid, '%e\t%e\t%e\t%e\t%e\t%e\t%e\t%s\n',results(i,1)/unitFactor,results(i,2)/unitFactor,results(i,3),results(i,4)*unitFactor,results(i,5)*unitFactor,results(i,6)*unitFactor,results(i,7),char(procfiles(i)));
% end
% fclose(fid);

end
