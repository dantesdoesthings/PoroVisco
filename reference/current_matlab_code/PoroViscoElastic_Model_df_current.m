% POROVISCOELASTIC MODEL FIT FOR NANOINDENTATION 
% Author : MRI (ECU, Oyen Lab) 07/30/2019
% Edited by DMF (Myers Lab, Columbia) Nov. 2021
% Edited by AMM (Oyen Lab, WashU) July 2022

% INPUT : % test data (DT) > format- time(s) load(N) disp(m)
% OUTPUT : PP > format > G0(kPa) E0 Ginf Einf Ginf/G0 Tau1(s) Tau2 RSq-visco G(kPa) E neu kappa(m^2) D (m^2/s) RSq-poro RSq-PVE Model (Combined) 
% OUTPUT : PD > format > DT fit_force
function [PP PD]=PoroViscoElastic_Model_df_AM(DT,R,rampTime,holdTime,plotflag,i,leg,dwellid)
               
%% RELAXATION REGIME DATA
hmax=abs(DT(dwellid,3));
Pmax=DT(dwellid,2);  

%% INDENTER PARAMETERS
% Spherical Indenter
m=3/2;                                    % h^exponent
R=R;                                      % indenter radius (m)    
av=8.*sqrt(R)/3.;
ap=sqrt(R*hmax);

%% LEAST SQUARE FIT - Hu et al. APL 96, 121904 (2010)

% Poroelastic Initialization
P0guess = DT(2,2)*1.03;                         %Multiple by factor for improved fitting of data
Pinfguess = DT(end,2)*0.95;
Gguess = 3*P0guess/(16*hmax*ap);
vguess = 1 - (P0guess/(2*Pinfguess));             %fixed from orig MRI code (previous vguess =0.2)
kappaGuess=1E-14;                        %kappa/eta
logKappaGuess =log10(kappaGuess);
Dguess = (2*(1-vguess)/(1-2*vguess))*Gguess*kappaGuess;
logDguess = log10(Dguess);

% VARIABLE BOUNDS- D P0 v
if vguess<0
vmin = 0.5/vguess;
vmax = -1/vguess;
else
vmax = 0.5/vguess;
vmin = -1/vguess;       %fixed from original MRI code 
end

LB1 = [0 0 vmin];     % lower bound 
UB1 = [1 1 vmax];     % upper bound 
X01=[1 1 1];
Xg1 = [logDguess, Pinfguess, vguess];

% Viscoelastic Initialization
numParam=2;
C0guess=DT(end,2)/(hmax^m*av);                    
Cguess = C0guess*ones(1,numParam);
Gguess=(C0guess+sum(Cguess))/2;
Tguess = [rampTime holdTime];
X02 = ones(1, 2 * numParam + 1);
Xg2 = [C0guess Cguess Tguess];
LB2 = [0 0 0 0 0];
UB2 = [1 1 1 1 1];


%% OPT RUN
X0=[X01 X02];
Xguess=[Xg1 Xg2];
LB=[LB1 LB2];
UB=[UB1 UB2];
options.optim = optimset('MaxFunEvals',500000,'Display','none','TolX',1E-30);
X = lsqnonlin(@OBJPoroVisco_mri,X0,LB,UB,options.optim,DT,Xguess,rampTime,hmax,av,m,R);

%% POROELASTIC PARAMETER ESTIMATION
logDid = real(X(1).*logDguess);
Did = 10^logDid;
Pinf = X(2)*Pinfguess;
v = X(3)*vguess;
P0=Pinf*(2*(1-v));
G = 3*P0/(16*hmax*ap); % unit (Pa)
kappa = Did*(1-2*v)/(2*(1-v)*G);
kappa=0.89*1E-3*kappa; 
E=2*G*(1+0.5);

%% VISCOELASTIC PARAMETER ESTIMATION
Xid = real(X(4:end).*Xguess(4:end));            
C0 = Xid(1);
C=Xid(2:numParam +1);
T=Xid(numParam + 2:2*numParam + 1);
RCF=T./rampTime.*(exp(rampTime./T)-1);
C1=C(1)/RCF(1); C2=C(2)/RCF(2);
G0=(C0+sum(C))/2; E0=2*G0*(1+0.5);          % nu=0.5, unit(Pa) 
Ginf=C0/2;        Einf=2*Ginf*(1+0.5);      % nu=0.5  unit(Pa)
Tau=T;


%% LOAD-RELAXATION DATA 
for k=length(DT(:,1)):-1:1
tau(k) = Did*DT(k,1)/(R*hmax);
g(k) = 0.491*exp(-0.908*sqrt(tau(k))) + 0.509*exp(-1.679*(tau(k)));
fitP(k) = g(k)*(P0 - Pinf) + Pinf; 
fitV(k)=C0+C1*exp(-DT(k,1)/T(1))+C2*exp(-DT(k,1)/T(2));
fitV(k)=fitV(k)*av*hmax^(3/2);
fit(k)=fitP(k)*fitV(k)/DT(end,2);
end

s=1:1:size(DT,1);
RSq_visco= cofDet(DT(:,2), fitV(s)');         % Goodness of Fit (R2) 
RSq_poro= cofDet(DT(:,2), fitP(s)');         % Goodness of Fit (R2)
RSq_both= cofDet(DT(:,2), fit(s)');         % Goodness of Fit (R2) 

PP1=[G0/1000 E0/1000 Ginf/1000 Einf/1000 Ginf/G0 Tau(1) Tau(2) RSq_visco];
PP2 =[G/1000 E/1000 v kappa Did RSq_poro RSq_both];
PP=[PP1 PP2];
PD=[DT fit'];

%% CHECK PLOT

if plotflag==1
figure(i);
plot(DT(s,1),DT(s,2),'ok','MarkerSize',5,'LineWidth',1); hold on
s=1:1:size(DT,1);
plot(DT(s,1),fit(s)','-r','LineWidth',1); hold on
plot(DT(s,1),fitV(s)','--b','LineWidth',1); hold on
plot(DT(s,1),fitP(s)','-.g','LineWidth',1); hold on
ylabel('Load, {\it P}','fontsize',12);
xlabel('Hold time, {\it t-t}_r (s)','fontsize',12);
title(leg(i), ' PoroVisco Model');
set(gcf,'color','white'); box on
%xlim([0 10])
%ylim([0 1]);
lh=legend('Exp','{\it P}_{PVE}','{\it P}_{VE}','{\it P}_{PE}');
DIMX = 3; DIMY =2.75;
set(gcf,'Units','inches');
set(gca,'Units','inches');
box on;
set(gca, 'FontName', 'Times New Roman','FontSize',12);
set(gca,'LineWidth',1);
hold off;
end

end

function errvect=OBJPoroVisco_mri(X,expdata,Xguess,riseTime,hmax,a,m,R)
    X1=X(1:3);          X2=X(4:end);
    Xg1=Xguess(1:3);    Xg2=Xguess(4:end);
    % Poroelastic
    logD = X1(1) * Xg1(1); 
    Pinf = X1(2)*Xg1(2);
    v = X1(3) * Xg1(3);
    D = 10^logD;
    P0=Pinf*(2*(1-v));
    % Viscoelastic
    numParam = (length(Xg2)-1)/2;
    C0 =X2(1)*Xg2(1);
    C = X2(2:numParam + 1).*Xg2(2:numParam + 1);
    T = X2(numParam + 2:length(Xg2)).*Xg2(numParam + 2:length(Xg2));
    RCF=T./riseTime.*(exp(riseTime./T)-1);
    C1=C(1)/RCF(1); C2=C(2)/RCF(2);

    for k=1:length(expdata(:,1))
        tau(k) = D*expdata(k,1)/(R*hmax);
        g(k) = 0.491*exp(-0.908*sqrt(tau(k))) + 0.509*exp(-1.679*(tau(k)));
        fitP(k) = g(k)*(P0 - Pinf) + Pinf;
        fitV(k)=C0+C1*exp(-expdata(k,1)/T(1))+C2*exp(-expdata(k,1)/T(2));
        fitV(k)=fitV(k)*a*hmax^(3/2); 
        fit(k)=fitV(k)*fitP(k)/expdata(end,2);
    end
    errvect=(fit'-expdata(:,2))/mean(expdata(:,2)); 
    
    G0=(C0+sum(C))/2;
    ap=sqrt(R*hmax);
    G = 3*P0/(16*hmax*ap); % unit (Pa)
    if ((G-G0)/G0)>0.05
    errvect=errvect*1E6;
    end
end

