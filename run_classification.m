function [S,fl,casestudy,nlf,Co,Tix,Tq,index] =  run_process(x,ncase,reply2)
%% Load Case study  
switch ncase
    case 1
        case_eightbus_Urdaneta_LP; casestudy=' Urdaneta/EzzaddineLP et al';
    case 2
        case_eightbus_Birla; casestudy=' Birla et al';
    case 3
        case_eightbus_Ezeddinne; casestudy=' Ezzeddine NLP et al';
    case 4
        case_eightbus_Mahari; casestudy=' Mahari et al';
    case 5
        case_eightbus_Alipour; casestudy='Alipour IGSO';
    case 6
        case_eightbus_Meskin; casestudy=' Meskin et al';
    case 7
        case_eightbus_Sorrentino_LPSI; casestudy=' Sorrentino_LPSI et al';
    case 8
        case_eightbus_Sorrentino_LPVI; casestudy=' Sorrentino_LPVI et al';
     case 9
        case_eightbus_Sorrentino_LPEI; casestudy=' Sorrentino_LPEI et al';
    case 10 
        case_eightbus_Oliveira_LP; casestudy='Oliveira et al';
    otherwise
        case_eightbus_Sorrentino_LPSI; casestudy=' Sorrentino_LPSI et al';
end

fl=zeros(16,1);flag=0;
for jj=1:nlf
L(1,jj)=x;% fault distance xmin < x < xmax
end  
[index]=run_shortcircuit(L(1,:),reply2);%Run external short-circuit program 
for k=1:length(index(:,1)) %number of relay pairs
beta(index(k,3),index(k,4))= index(k,13);  %beta i main first interval due to fault in line k 
beta(index(k,1),index(k,2))= index(k,14);  %beta j backup first interval due to fault in line k   
betap(index(k,3),index(k,4))= index(k,15);  %beta i' main transient due to fault in line k when q is open
betap(index(k,1),index(k,2))= index(k,16);  %beta j' backup transient due to fault in line k when q is open
theta(index(k,3),index(k,4))= index(k,17);  %beta i main first interval due to fault in line k 
theta(index(k,1),index(k,2))= index(k,18);  %beta j backup first interval due to fault in line k   
thetap(index(k,3),index(k,4))= index(k,19);  %beta i' main transient due to fault in line k when q is open
thetap(index(k,1),index(k,2))= index(k,20);  %beta j' backup transient due to fault in line k when q is open
gammap(index(k,22),index(k,1),index(k,3),index(k,2))= index(k,21);% gamma' q j i k
gammapp(index(k,22),index(k,1),index(k,3),index(k,2))= index(k,24);% gamma'' q j i k  relay j
gammappp(index(k,22),index(k,1),index(k,3),index(k,2))= index(k,25);% gamma''' q j i k relaj i
end 
i=zeros(1,length(index(:,1)));
q=zeros(1,length(index(:,1)));
p=zeros(2,length(index(:,1)));
j=zeros(2,length(index(:,1)));
Tix=zeros(1,length(index(:,1)));
Tq=zeros(1,length(index(:,1)));
for k=1:length(index(:,1))
 if   beta(index(k,22),index(k,4))*D(index(k,22)) > 0
 if  beta(index(k,3),index(k,4))*D(index(k,3)) > beta(index(k,22),index(k,4))*D(index(k,22)) %Who trip first
 q(k)=index(k,22);
 i(k)=index(k,3);
 J=find(index(:,22)==q(k));
 for ii=1:length(J)
     j(ii,k)=index(J(ii),1);
 end
  P=find(index(:,3)==q(k));
 for ii=1:length(P)
     p(ii,k)=index(P(ii),1);
 end
 end
 else
 q(k)=index(k,3);
 i(k)=index(k,22);
 J=find(index(:,3)==q(k));
 for ii=1:length(J)
     j(ii,k)=index(J(ii),1);
 end
  P=find(index(:,22)==q(k));
 for ii=1:length(P)
     p(ii,k)=index(P(ii),1);
 end
end
end % Identifies for every fault who is q, i, j and p
Ad=vertcat(i,j);
Bd=unique(Ad','rows')';
Bd(:,1) = [];
Ai=vertcat(q,p);
Bi=unique(Ai','rows')';
Bi(:,1) = [];
iter=1;
for h=1:nlf
for m=2:3
if  Bd(m,h) > 0 
C1(iter,1)=Bd(m,h);
C1(iter,2)=Bd(1,h);
C1(iter,3)=h;
C1(iter,4)=Bi(1,h);
iter=iter+1;
end
end
end%C1 indicates tripping sequence ij 
iter=1;
for h=1:nlf
for m=2:3
if  Bi(m,h) > 0 
C2(iter,1)=Bi(m,h);
C2(iter,2)=Bi(1,h);
C2(iter,3)=h;
iter=iter+1;
end
end
end%C2 indicates tripping sequence qp
%% Separation times                        
for g=1:length(C2(:,1))
if abs(theta(C2(g,1),C2(g,3))-theta(C2(g,2),C2(g,3))) < qmax &...%period 1: no reverse current relay p 
        beta(C2(g,1),C2(g,3)) > 0 %period 1: no loss of sensitivity relay p  
flag=flag+1; 
fl(1)=fl(1)+1; 
S(flag)=beta(C2(g,1),C2(g,3))*D(C2(g,1))-beta(C2(g,2),C2(g,3))*D(C2(g,2));
Tq(flag)=beta(C2(g,2),C2(g,3))*D(C2(g,2));
end
end% Type 1  - Case 1  
for g=1:length(C1(:,1))
if abs(theta(C1(g,1),C1(g,3))-theta(C1(g,2),C1(g,3))) < qmax &... %period 1: no reverse current relay j
   abs(thetap(C1(g,1),C1(g,3))-thetap(C1(g,2),C1(g,3))) < qmax &... %period 2: no reverse current relay j
   beta(C1(g,1),C1(g,3)) > 0 &...%period 1: no loss of sensitivity relay j
   betap(C1(g,1),C1(g,3)) > 0 &...%period 2: no loss of sensitivity relay j
   beta(C1(g,2),C1(g,3)) > 0 &...%period 1: no loss of sensitivity relay i
   betap(C1(g,2),C1(g,3)) > 0  %period 2: no loss of sensitivity relay i
flag=flag+1; 
fl(2)=fl(2)+1;
S(flag)=betap(C1(g,1),C1(g,3))*D(C1(g,1))-betap(C1(g,2),C1(g,3))*D(C1(g,2))-gammap(C1(g,4),C1(g,1),C1(g,2),C1(g,3))*D(C1(g,4));
% betap(C1(g,1),C1(g,3)) 
% betap(C1(g,2),C1(g,3))
% gammap(C1(g,4),C1(g,1),C1(g,2),C1(g,3))
% D(C1(g,1))
% D(C1(g,2))
% D(C1(g,4))
% % C1(g,1)
% % C1(g,3)
% % C1(g,2) 
% % C1(g,3)
% % % C1(g,4) 
% % % C1(g,1)
% % % C1(g,2)
% % % C1(g,3)
% % % 
%  pause
Tj=beta(C1(g,1),C1(g,3))*D(C1(g,1));
Ti=beta(C1(g,2),C1(g,3))*D(C1(g,2));
Tqq=beta(C1(g,4),C1(g,3))*D(C1(g,4));
Tpj=betap(C1(g,1),C1(g,3))*D(C1(g,1));
Tpi=betap(C1(g,2),C1(g,3))*D(C1(g,2));
Tix(flag)=Tqq+Tpi*(1-Tqq/Ti);
Tjx(flag)=Tqq+Tpj*(1-Tqq/Tj);
end
end% Type 2  - Case 2   
for g=1:length(C1(:,1))
if abs(theta(C1(g,1),C1(g,3))-theta(C1(g,2),C1(g,3))) < qmax &... %period 1: no reverse current relay j
   abs(thetap(C1(g,1),C1(g,3))-thetap(C1(g,2),C1(g,3))) < qmax &... %period 2: no reverse current relay j
   beta(C1(g,1),C1(g,3)) < 0 &...%period 1: yes loss of sensitivity relay j
   betap(C1(g,1),C1(g,3)) > 0 &...%period 2: no loss of sensitivity relay j
   beta(C1(g,2),C1(g,3)) > 0 &...%period 1: no loss of sensitivity relay i
   betap(C1(g,2),C1(g,3)) > 0%period 2: no loss of sensitivity relay i

flag=flag+1;
fl(3)=fl(3)+1; 
S(flag)=betap(C1(g,1),C1(g,3))*D(C1(g,1))-betap(C1(g,2),C1(g,3))*D(C1(g,2))-gammapp(C1(g,4),C1(g,1),C1(g,2),C1(g,3))*D(C1(g,4));
Tj=beta(C1(g,1),C1(g,3))*D(C1(g,1));
Ti=beta(C1(g,2),C1(g,3))*D(C1(g,2));
Tqq=beta(C1(g,4),C1(g,3))*D(C1(g,4));
Tpj=betap(C1(g,1),C1(g,3))*D(C1(g,1));
Tpi=betap(C1(g,2),C1(g,3))*D(C1(g,2));
Tix(flag)=Tqq+Tpi*(1-Tqq/Ti);
Tjx(flag)=Tix(flag)+S(flag);
end
end% Type 3a - Case 3   
for g=1:length(C1(:,1))
if abs(theta(C1(g,1),C1(g,3))-theta(C1(g,2),C1(g,3))) > qmax &... %period 1: yes reverse current relay j
   abs(thetap(C1(g,1),C1(g,3))-thetap(C1(g,2),C1(g,3))) < qmax &... %period 2: no reverse current relay j
   beta(C1(g,1),C1(g,3)) > 0 &...%period 1: no loss of sensitivity relay j
   betap(C1(g,1),C1(g,3)) > 0 &...%period 2: no loss of sensitivity relay j
   beta(C1(g,2),C1(g,3)) > 0 &...%period 1: no loss of sensitivity relay i
   betap(C1(g,2),C1(g,3)) > 0 %period 2: no loss of sensitivity relay i
    
flag=flag+1; 
fl(4)=fl(4)+1; 
S(flag)=betap(C1(g,1),C1(g,3))*D(C1(g,1))-betap(C1(g,2),C1(g,3))*D(C1(g,2))-gammapp(C1(g,4),C1(g,1),C1(g,2),C1(g,3))*D(C1(g,4));
Tj=beta(C1(g,1),C1(g,3))*D(C1(g,1));
Ti=beta(C1(g,2),C1(g,3))*D(C1(g,2));
Tqq=beta(C1(g,4),C1(g,3))*D(C1(g,4));
Tpj=betap(C1(g,1),C1(g,3))*D(C1(g,1));
Tpi=betap(C1(g,2),C1(g,3))*D(C1(g,2));
Tix(flag)=Tqq+Tpi*(1-Tqq/Ti);
Tjx(flag)=Tix(flag)+S(flag);
 
end
end% Type 3b - Case 4   
for g=1:length(C1(:,1))
if abs(theta(C1(g,1),C1(g,3))-theta(C1(g,2),C1(g,3))) > qmax &... %period 1: yes reverse current relay j
   abs(thetap(C1(g,1),C1(g,3))-thetap(C1(g,2),C1(g,3))) < qmax &... %period 2: no reverse current relay j
   beta(C1(g,1),C1(g,3)) < 0 &...%period 1: yes loss of sensitivity relay j
   betap(C1(g,1),C1(g,3)) > 0 &...%period 2: no loss of sensitivity relay j
   beta(C1(g,2),C1(g,3)) > 0 &...%period 1: no loss of sensitivity relay i
   betap(C1(g,2),C1(g,3)) > 0 %period 2: no loss of sensitivity relay iflag=flag+1; 
fl(5)=fl(5)+1; 
S(flag)=betap(C1(g,1),C1(g,3))*D(C1(g,1))-betap(C1(g,2),C1(g,3))*D(C1(g,2))-gammapp(C1(g,4),C1(g,1),C1(g,2),C1(g,3))*D(C1(g,4));
Tj=beta(C1(g,1),C1(g,3))*D(C1(g,1));
Ti=beta(C1(g,2),C1(g,3))*D(C1(g,2));
Tqq=beta(C1(g,4),C1(g,3))*D(C1(g,4));
Tpj=betap(C1(g,1),C1(g,3))*D(C1(g,1));
Tpi=betap(C1(g,2),C1(g,3))*D(C1(g,2));
Tix(flag)=Tqq+Tpi*(1-Tqq/Ti);
Tjx(flag)=Tix(flag)+S(flag);
end
end% Type 3c - Case 5  
for g=1:length(C1(:,1))
if abs(theta(C1(g,1),C1(g,3))-theta(C1(g,2),C1(g,3))) < qmax &... %period 1: no reverse current relay j
   abs(thetap(C1(g,1),C1(g,3))-thetap(C1(g,2),C1(g,3))) < qmax &... %period 2: no reverse current relay j
   beta(C1(g,1),C1(g,3)) > 0 &...%period 1: no loss of sensitivity relay j
   betap(C1(g,1),C1(g,3)) > 0 &...%period 2: no loss of sensitivity relay j
   beta(C1(g,2),C1(g,3)) < 0 &...%period 1: yes loss of sensitivity relay i
   betap(C1(g,2),C1(g,3)) >0  %period 2: no loss of sensitivity relay iif abs(theta(C1(g,1),C1(g,3))-theta(C1(g,2),C1(g,3))) < qmax &  abs(thetap(C1(g,1),C1(g,3))-thetap(C1(g,2),C1(g,3))) < qmax & beta(C1(g,1),C1(g,3)) > 0 & betap(C1(g,1),C1(g,3)) > 0 & beta(C1(g,2),C1(g,3)) < 0 & betap(C1(g,2),C1(g,3)) > 0  & beta(C1(g,4),C1(g,3))*D(C1(g,4)) 
flag=flag+1; 
fl(6)=fl(6)+1; 
S(flag)=betap(C1(g,1),C1(g,3))*D(C1(g,1))-betap(C1(g,2),C1(g,3))*D(C1(g,2))-gammappp(C1(g,4),C1(g,1),C1(g,2),C1(g,3))*D(C1(g,4));
Tj=beta(C1(g,1),C1(g,3))*D(C1(g,1));
Ti=beta(C1(g,2),C1(g,3))*D(C1(g,2));
Tqq=beta(C1(g,4),C1(g,3))*D(C1(g,4));
Tpj=betap(C1(g,1),C1(g,3))*D(C1(g,1));
Tpi=betap(C1(g,2),C1(g,3))*D(C1(g,2));
Tjx(flag)=Tqq+Tpj*(1-Tqq/Tj);
Tix(flag)=Tjx(flag)-S(flag);
end
end% Type 4  - Case 6  
for g=1:length(C1(:,1))
if abs(theta(C1(g,1),C1(g,3))-theta(C1(g,2),C1(g,3))) < qmax &... %period 1: no reverse current relay j
   abs(thetap(C1(g,1),C1(g,3))-thetap(C1(g,2),C1(g,3))) < qmax &... %period 2: no reverse current relay j
   beta(C1(g,1),C1(g,3)) < 0 &...%period 1: yes loss of sensitivity relay j
   betap(C1(g,1),C1(g,3)) > 0 &...%period 2: no loss of sensitivity relay j
   beta(C1(g,2),C1(g,3)) < 0 &...%period 1: yes loss of sensitivity relay i
   betap(C1(g,2),C1(g,3)) > 0  %period 2: no loss of sensitivity relay i
fl(7)=fl(7)+1; 
S(flag)=betap(C1(g,1),C1(g,3))*D(C1(g,1))-betap(C1(g,2),C1(g,3))*D(C1(g,2))-0*gammappp(C1(g,4),C1(g,1),C1(g,2),C1(g,3))*D(C1(g,4));
Tj=beta(C1(g,1),C1(g,3))*D(C1(g,1));
Ti=beta(C1(g,2),C1(g,3))*D(C1(g,2));
Tqq=beta(C1(g,4),C1(g,3))*D(C1(g,4));
Tpj=betap(C1(g,1),C1(g,3))*D(C1(g,1));
Tpi=betap(C1(g,2),C1(g,3))*D(C1(g,2));
Tix(flag)=Tpi;
Tjx(flag)=Tpj;
end
end% Type 5a - Case 7   
for g=1:length(C1(:,1))
if abs(theta(C1(g,1),C1(g,3))-theta(C1(g,2),C1(g,3))) > qmax &... %period 1: yes reverse current relay j
   abs(thetap(C1(g,1),C1(g,3))-thetap(C1(g,2),C1(g,3))) < qmax &... %period 2: no reverse current relay j
   beta(C1(g,1),C1(g,3)) > 0 &...%period 1: no loss of sensitivity relay j
   betap(C1(g,1),C1(g,3)) > 0 &...%period 2: no loss of sensitivity relay j
   beta(C1(g,2),C1(g,3)) < 0 &...%period 1: yrs loss of sensitivity relay i
   betap(C1(g,2),C1(g,3)) >0  %period 2: no loss of sensitivity relay iif abs(theta(C1(g,1),C1(g,3))-theta(C1(g,2),C1(g,3))) > qmax &  abs(thetap(C1(g,1),C1(g,3))-thetap(C1(g,2),C1(g,3))) < qmax & beta(C1(g,1),C1(g,3)) > 0 & betap(C1(g,1),C1(g,3)) > 0 & beta(C1(g,2),C1(g,3)) < 0 & betap(C1(g,2),C1(g,3)) > 0  & beta(C1(g,4),C1(g,3))*D(C1(g,4)) 
flag=flag+1; 
fl(8)=fl(8)+1; 
S(flag)=betap(C1(g,1),C1(g,3))*D(C1(g,1))-betap(C1(g,2),C1(g,3))*D(C1(g,2))-0*gammappp(C1(g,4),C1(g,1),C1(g,2),C1(g,3))*D(C1(g,4));
Tj=beta(C1(g,1),C1(g,3))*D(C1(g,1));
Ti=beta(C1(g,2),C1(g,3))*D(C1(g,2));
Tqq=beta(C1(g,4),C1(g,3))*D(C1(g,4));
Tpj=betap(C1(g,1),C1(g,3))*D(C1(g,1));
Tpi=betap(C1(g,2),C1(g,3))*D(C1(g,2));
Tix(flag)=Tpi;
Tjx(flag)=Tpj;
end
end% Type 5b - Case 8   
for g=1:length(C1(:,1))
if abs(theta(C1(g,1),C1(g,3))-theta(C1(g,2),C1(g,3))) > qmax &... %period 1: yes reverse current relay j
   abs(thetap(C1(g,1),C1(g,3))-thetap(C1(g,2),C1(g,3))) < qmax &... %period 2: no reverse current relay j
   beta(C1(g,1),C1(g,3)) < 0 &...%period 1: yes loss of sensitivity relay j
   betap(C1(g,1),C1(g,3)) > 0 &...%period 2: no loss of sensitivity relay j
   beta(C1(g,2),C1(g,3)) < 0 &...%period 1: yes loss of sensitivity relay i
   betap(C1(g,2),C1(g,3)) >0  %period 2: no loss of sensitivity relay iif abs(theta(C1(g,1),C1(g,3))-theta(C1(g,2),C1(g,3))) > qmax &  abs(thetap(C1(g,1),C1(g,3))-thetap(C1(g,2),C1(g,3))) < qmax & beta(C1(g,1),C1(g,3)) > 0 & betap(C1(g,1),C1(g,3)) > 0 & beta(C1(g,2),C1(g,3)) < 0 & betap(C1(g,2),C1(g,3)) > 0  & beta(C1(g,4),C1(g,3))*D(C1(g,4)) 
flag=flag+1; 
fl(9)=fl(9)+1; 
S(flag)=betap(C1(g,1),C1(g,3))*D(C1(g,1))-betap(C1(g,2),C1(g,3))*D(C1(g,2))-0*gammappp(C1(g,4),C1(g,1),C1(g,2),C1(g,3))*D(C1(g,4));
Tj=beta(C1(g,1),C1(g,3))*D(C1(g,1));
Ti=beta(C1(g,2),C1(g,3))*D(C1(g,2));
Tqq=beta(C1(g,4),C1(g,3))*D(C1(g,4));
Tpj=betap(C1(g,1),C1(g,3))*D(C1(g,1));
Tpi=betap(C1(g,2),C1(g,3))*D(C1(g,2));
Tix(flag)=Tpi;
Tjx(flag)=Tpj;
end
end% Type 5c - Case 9   
for g=1:length(C2(:,1))
if abs(theta(C2(g,1),C2(g,3))-theta(C2(g,2),C2(g,3))) > qmax &...period 1:yes reverse current relay p 
   beta(C2(g,1),C2(g,3)) > 0   %period 1: no loss of sensitivity relay p  
fl(10)=fl(10)+1; 
flag=flag+1; 
Tq(flag)=beta(C2(g,2),C2(g,3))*D(C2(g,2));
end
end% Type 6a - Case 10  
for g=1:length(C1(:,1))
if abs(theta(C1(g,1),C1(g,3))-theta(C1(g,2),C1(g,3))) > qmax &... %period 1: yes reverse current relay j
   abs(thetap(C1(g,1),C1(g,3))-thetap(C1(g,2),C1(g,3))) > qmax &... %period 2:yes reverse current relay j
   beta(C1(g,1),C1(g,3)) > 0 &...%period 1: no loss of sensitivity relay j
   betap(C1(g,1),C1(g,3)) > 0 &...%period 2: no loss of sensitivity relay j
   beta(C1(g,2),C1(g,3)) > 0 &...%period 1: no loss of sensitivity relay i
   betap(C1(g,2),C1(g,3)) > 0 %period 2: no loss of sensitivity relay i
fl(11)=fl(11)+1; 
end
end% Type 6b - Case 11   
for g=1:length(C2(:,1))
if abs(theta(C2(g,1),C2(g,3))-theta(C2(g,2),C2(g,3))) < qmax &...period 1:no reverse current relay p 
        beta(C2(g,1),C2(g,3)) < 0  %period 1: yes  loss of sensitivity relay p
fl(12)=fl(12)+1;
flag=flag+1; 
Tq(flag)=beta(C2(g,2),C2(g,3))*D(C2(g,2));
end
end% Type 6c - Case 12  
for g=1:length(C2(:,1))
if abs(theta(C2(g,1),C2(g,3))-theta(C2(g,2),C2(g,3))) > qmax &...period 1:yes reverse current relay p 
        beta(C2(g,1),C2(g,3)) < 0  %period 1: yes  loss of sensitivity relay p
fl(13)=fl(13)+1; 
flag=flag+1; 
Tq(flag)=beta(C2(g,2),C2(g,3))*D(C2(g,2));
end
end% Type 6d - Case 13  
for g=1:length(C1(:,1))
if    beta(C1(g,2),C1(g,3)) < 0 &...%period 1: yes loss of sensitivity relay i
   betap(C1(g,2),C1(g,3)) < 0 %period 2: yes loss of sensitivity relay i
 fl(14)=fl(14)+1; 
 end
end% Type 6e - Case 14   
for g=1:length(C1(:,1))
if  betap(C1(g,1),C1(g,3)) < 0  %period 2: yes loss of sensitivity relay j
fl(15)=fl(15)+1; 
end
end% Type 6f - Case 15  
for g=1:length(C1(:,1))
if abs(theta(C1(g,1),C1(g,3))-theta(C1(g,2),C1(g,3))) > qmax &... %period 1: yes reverse current relay j
   abs(thetap(C1(g,1),C1(g,3))-thetap(C1(g,2),C1(g,3))) > qmax %period 2: yes reverse current relay j    
fl(16)=fl(16)+1; 
end
end% Type 6g - Case 16  
 



