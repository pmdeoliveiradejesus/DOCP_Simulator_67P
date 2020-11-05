% DOCR (67 PHASE) - Directional 67 PHASE Overcurrent Coordination Relays Problem
% This program simulates the response of the protection system (speed & selectivity) 
% considering transient configurations
% Developer:
% Paulo M. De Oliveira, pm.deoliveiradejes@uniandes.edu.co
% https://power.uniandes.edu.co/
% Los Andes University
% First version: November 6, 2020 
% Documentation:
% De Oliveira-De Jesus, P.M. and Sorrentino E. 
% Methodology to assess global performance indexes for sensitivity, selectivity, and speed of directional overcurrent protection systems
% submitted to Electrical Power Systems Research

clear all
close all
clc
time000=cputime;
%% Initial screen
% disp('*******************************************************')
 disp('DOCR (67 Phase)                                        ')
 disp('67 Relay System Simulator                              ')
 disp('Version 1.0 (c) 2020                                   ')
 disp('Power and Energy Group - https://power.uniandes.edu.co/')
 disp('Universidad de los Andes, Colombia')    
 disp('*******************************************************')
addpath('./data/')
%% Available case studies
% 1: Ezzeddine/Urdaneta LP	
% 2: Birla	
% 3: Ezzeddine NLP	
% 4: Mahari	
% 5: Alipour
% 6: Meskin	
% 7: Sorrentino
% 8: De Oliveira
ncase=7;% Case number selected
disp('Available study cases:')
disp('1: Ezzeddine/Urdaneta')
disp('2: Birla')
disp('3: Ezzeddine')
disp('4: Mahari')
disp('5: Alipour')
disp('6: Meskin')
disp('7: Sorrentino') 
reply = input('Please select the case number from the list 1-8 [7]: ');
if isempty(reply)
          ncase = 7;
else
          ncase=reply;
end 
reply2 = input('Do you want to include a prefault power flow? y/n [n]:','s');
if isempty(reply2)
reply2 = 'n';
end         
neval=1000;% Number of faults per line uniformily distributed
Case0=zeros(16,1);% Initialize type pairs vector
T=[0];%Initialize primary times
i=1;%Initialize flag for time progress
iter=1;%initialize counter for separation times
%% Begins Iterative process 
for k=1:neval  
i=i+1;
if i==neval/4
fprintf('Wait! simulating... Progress: 25%%')
elseif i==neval/2
fprintf(', 50%%')
elseif i==3*neval/4
fprintf(', 75%%')    
elseif i==neval
fprintf(', 100%%\n')
end 
%x= 0.5; %with neval=1  x can be located in any place 0<x<1 (Optional)
x=k/(neval+1);%Uniform distributed faults, if neval=1000 x goes from 0.001 to 0.999 
[S,Case,casestudy,nlf,Co,Tix,Tq,index]=run_process(x,ncase,reply2);%Runs the classification script

Case0=Case0+Case;%All 16 pair types classified are aggregated here
T=unique(vertcat(unique([Tix';Tq']),T));%All nr primary times are aggregated here
sizeS(k)=length(S);%Set length of each separation vector
for k=1:length(S)
SepTime(iter)=S(k);
   iter=iter+1;
    end%All calculated separation times are aggregated here
end 
%% Iterative process ends
elapsedtime000=cputime-time000;% Set simulation time
ki=0;
for k=1:length(SepTime)
if SepTime(k) < Co
ki=ki+1;       
end
end%determine number of separation times below specified Co (CTI)
% Types 1 to 6 calculation
result(1,1)= Case0(1); %Number of relay pairs Type 1
result(2,1)= Case0(2); %Number of relay pairs Type 2
result(3,1)= Case0(3)+Case0(4)+Case0(5); %Number of relay pairs Type 3
result(4,1)= Case0(6); %Number of relay pairs Type 4
result(5,1)= Case0(7)+Case0(8)+Case0(9); %Number of relay pairs Type 5
result(6,1)= Case0(10)+Case0(11)+Case0(12)+Case0(13)+Case0(14)+Case0(15)+Case0(16); %Number of relay pairs Type 6 
Nf=result(1,1)+result(2,1)+result(3,1)+result(4,1)+result(5,1);% Number of calculable sep time backup-main relay pair  
Nnf=result(6,1);% Number of Non-Feasible relay pairs 
N=Nf+Nnf;% Total pairs
Nrev=Case0(4)+Case0(10)+Case0(11); % Pairs with reverse fault currents
Nnosen=N-Nrev-Case0(1)-Case0(2);% Pairs with loss of sensitivity
Nnosel=ki;% Pairs with loss of selectivity
%% Performance indexes
% selectivity
sel=(1-Nnosel/Nf)*100;%selectivity level index
minSepTime=(min(SepTime));%Minimum separation time (seconds)
meanSepTime=mean(SepTime);%Mean Separation Time (seconds)
devSepTime=std(SepTime);%Std. Dev. Separation Time (seconds)
% sensitivity
sen=100*(1-Nnosen/(N-Nrev));%sensitivity level index
% speed
T(T==0) = [];
meanPrimTime=mean(T);% Average primary operation time (seconds)
AvgPrimSpeed=1/meanPrimTime; %Average primary speed (1/seconds)
devPrimTime=std(T);% Average primary operation time (seconds)
meanBackTime=meanPrimTime+meanSepTime;% Average backup operation time (seconds)
%% Screen output
disp('*******************************************************')
fprintf('Case study:%s\n',casestudy)
fprintf('Simulation results:\n') 
fprintf('Prefault power flow included?: %s\n',reply2) 
fprintf('Number of fauls per line %d\n',neval)
fprintf('Simulated primary-backup pairs %d\n',(N))
fprintf('___________________________________________________________________________________\n');
fprintf('Relay response classification:\n');
fprintf('Type 1 Normal operation pairs j-i  :  %4d,  %4.1f %%\n',result(1,1), 100*result(1,1)/(N));
fprintf('Type 2 Normal operation pairs p-q  :  %4d,  %4.1f %%\n',result(2,1), 100*result(2,1)/(N));
fprintf('Type 3 Partial operation relay  j  :  %4d,  %4.1f %%\n',result(3,1), 100*result(3,1)/(N));
fprintf('Type 4 Partial operation relay  i  :  %4d,  %4.1f %%\n',result(4,1), 100*result(4,1)/(N));
fprintf('Type 5 Partial operation relays j-i:  %4d,  %4.1f %%\n',result(5,1), 100*result(5,1)/(N));
fprintf('Type 6 No operation                :  %4d,  %4.1f %%\n',result(6,1), 100*result(6,1)/(N));
fprintf('___________________________________________________________________________________\n');
fprintf('Feasible back-main pairs (there is a calculated separation time) :  %4d,  %4.1f %%\n',Nf, 100*Nf/(N));
fprintf('Non Feasible back-main pairs (there is no separation time S)     :  %4d,  %4.1f %%\n',Nnf, 100*Nnf/(N));
fprintf('Relay pairs without loss of sensitivity                          :  %4d,  %4.1f %%\n',N-Nrev, sen);
fprintf('Mean   Operation Times calc with specified TDSs and Ips          : %7.4f (ms)\n',meanPrimTime*1000)
fprintf('StdDev Operation Times calc with specified TDSs and IPs          : %7.4f (ms)\n',devPrimTime*1000)
fprintf('Verifying if minimum calculated separation time %5.4f is higher than the coordination time C=%7.4f\n',minSepTime,Co) 
fprintf('Result: %5.4f %% of the faults accomplishes the allowable coordination interval C=%7.4f\n',sel,Co) 
fprintf('Elapsed simulation time: %6.2f s \n',elapsedtime000)
disp('****************************************************************************************')
close all
% Figure 1 - Histogram of Separation Times
figure('name','Separation Times Histogram (s)','position',[0, 300, 400, 200])
xbins = -0.50:0.01:2;
h=histogram(SepTime,xbins,'FaceColor','yellow');
set(gcf,'color','w')
counts = h.Values';
% Figure 1 - Histogram of System Primary Relay Operation Times
figure('name','SPROT Histogram (s)','position',[0, 0, 400, 200])
% SROT Histogram
h=histogram(T,200,'FaceColor','yellow');
set(gcf,'color','w')
counts = h.Values;
% Save results for latex table
res(1,1)=N; % Total pairs
res(2,1)=Nf; % Total feasible pairs
res(3,1)=Nnf; % Total non-feasible pairs
res(4,1)=sen; % percentage
res(5,1)=sel; % percentage
res(6,1)=AvgPrimSpeed;%1/s
res(7,1)=meanPrimTime*1000; %ms
res(8,1)=devPrimTime*1000;%ms
res(9,1)=meanSepTime*1000;%ms
res(10,1)=minSepTime*1000;%ms
res(11,1)=elapsedtime000;% CPU time
Case0=Case0';
