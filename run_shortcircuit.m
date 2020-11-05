%Short-circuit program - Stevenson's Zbus apprroach
function [index] =  run_shortcircuit(L,reply2)
global bdat ldat tdat k1 k2 k3 Ip Ibase econv itermax 
% Relaying data
 nbuses=length(bdat(:,1));
    flag1=0;
    for k=1:length(bdat(:,1))
    if bdat(k,3) == 3
flag1=flag1+1;
    end
    end
ngen=length(bdat(:,1))-flag1;
nlines=length(ldat(:,1))-ngen;
bdat1=bdat;   
ldat1=ldat; 
for k=ngen+1:nbuses
    bdat1(k,:)=bdat(k,:);    
end
for k=nbuses+1:nlines+nbuses
    bdat1(k,1)=k;     
    bdat1(k,2)=999;    
    bdat1(k,3)=3;  
end    
for k=ngen+1:2*nlines+ngen
    ldat1(k,3)=k-ngen;        
end     
for k=1:nlines
ldat1(k*2+(ngen-1),1)=ldat(k+ngen,1);
end
for k=1:nlines
ldat1(1+k*2+(ngen-1),2)=ldat(k+ngen,2);
end
for k=1:nlines
ldat1(k*2+(ngen),1)=k+ngen+flag1;
ldat1(k*2+(ngen),4)=ldat(ngen+k,4)*(1-L(k));   
ldat1(k*2+(ngen),5)=ldat(ngen+k,5)*(1-L(k)); 
end
for k=1:nlines
ldat1(k*2+(ngen-1),2)=k+ngen+flag1;
ldat1(k*2+(ngen-1),4)=ldat(ngen+k,4)*(L(k));   
ldat1(k*2+(ngen-1),5)=ldat(ngen+k,5)*(L(k)); 
end
%Relay table
for i=1:length(ldat1(:,1))-ngen
    Relay(i,1)=i;%relay number
    Relay(i,2)=ldat1(i+ngen,1)-ngen; %from bus ns 
    Relay(i,3)=ldat1(i+ngen,2)-ngen; % to bus nr 
end
% Ordering relay table
for k=1:length(Relay(:,1))
if Relay(k,2) >  Relay(k,3)
flag5=Relay(k,2);
  Relay(k,2)=Relay(k,3);
  Relay(k,3)=flag5;
end
end
for k=ngen+1:length(ldat1(:,1))
ldat2(k-ngen,:)=ldat1(k,:);
ldat2(k-ngen,1)=ldat1(k,1)-ngen;
ldat2(k-ngen,2)=ldat1(k,2)-ngen;
end
%Index=zeros(10,30);
w=ones(30,30);
for  kk=flag1+1:length(bdat1(:,1))-ngen
for j=1:length(ldat2(:,1))
if ldat2(j,1)==kk
    extreme1(w(kk,1),1)=ldat2(j,2);
    w(kk,1)=w(kk,1)+1;
end
if ldat2(j,2)==kk
    extreme2(w(kk,2),1)=ldat2(j,1);
    w(kk,2)=w(kk,2)+1;
end
end
for jj=1:length(extreme1(:,1))
for j=1:length(ldat2(:,1))
if ldat2(j,1)==extreme1(jj,1) & ldat2(j,2)~=kk
    extreme1a(w(kk,3),1)=ldat2(j,2);
    w(kk,3)=w(kk,3)+1;
end
if ldat2(j,2)==extreme1(jj,1) & ldat2(j,1)~=kk
    extreme1a(w(kk,3),1)=ldat2(j,1);
    w(kk,3)=w(kk,3)+1;
end
end
end
for jj=1:length(extreme1a(:,1))
for j=1:length(ldat2(:,1))
if ldat2(j,1)==extreme1a(jj,1) %& ldat2(j,2)~=extreme1(1,1)
    extreme1b(w(kk,4),1)=ldat2(j,2);
    w(kk,4)=w(kk,4)+1;
end
if ldat2(j,2)==extreme1a(jj,1) %& ldat2(j,1)~=extreme1(1,1)
    extreme1b(w(kk,4),1)=ldat2(j,1);
    w(kk,4)=w(kk,4)+1;
end
end
end
for k=1:length(extreme1b)
    if extreme1b(k,1) == extreme1(1,1)
    else
       extreme1e(w(kk,5),1)=extreme1b(k,1);
       w(kk,5)= w(kk,5)+1;
    end
end
 
for jj=1:length(extreme2(:,1))
for j=1:length(ldat2(:,1))
if ldat2(j,1)==extreme2(jj,1) & ldat2(j,2)~=kk
    extreme2a(w(kk,6),1)=ldat2(j,2);
    w(kk,6)=w(kk,6)+1;
end
if ldat2(j,2)==extreme2(jj,1) & ldat2(j,1)~=kk
    extreme2a(w(kk,6),1)=ldat2(j,1);
    w(kk,6)=w(kk,6)+1;
end
end
end
for jj=1:length(extreme2a(:,1))
for j=1:length(ldat2(:,1))
if ldat2(j,1)==extreme2a(jj,1) 
    extreme1c(w(kk,7),1)=ldat2(j,2);
    w(kk,7)=w(kk,7)+1;
end
if ldat2(j,2)==extreme2a(jj,1) 
    extreme1c(w(kk,7),1)=ldat2(j,1);
    w(kk,7)=w(kk,7)+1;
end
end
end
for k=1:length(extreme1c)
    if extreme1c(k,1) == extreme2(1,1)
    else
       extreme1d(w(kk,8),1)=extreme1c(k,1);
       w(kk,8)=w(kk,8)+1;
    end
end
InM=[extreme2,extreme1];
InB2=[length(extreme1e),extreme1e',extreme1a']; 
InB1=[length(extreme1d),extreme1d',extreme2a']; 
for k=1:length(InM)
    IM(kk,k)=InM(k);  
end
for k=1:length(InB2)
    IB2(kk,k)=InB2(k);  
end
for k=1:length(InB1)
    IB1(kk,k)=InB1(k);  
end
clear extreme1 extreme2 extreme1b extreme1c extreme1e extreme1a extreme1d extreme2a
end
 cc=length(bdat1(:,1))-ngen;
 dd=length(ldat2(:,1));
  for k=1:ngen
ldat2(dd+k,:)=[cc+k,ldat(k,2)-ngen,11113,tdat(k,4),tdat(k,5),0,0,0];
  end
 aa=length(bdat1(:,1));
 bb=length(ldat1(:,1));
 for k=1:ngen
 ldat1(k,2)=aa+k;
 end
  for k=1:ngen
  ldat1(bb+k,:)=[aa+k,ldat(k,2),111113,tdat(k,4),tdat(k,5),0,0,0];
  bdat1(aa+k,:)=[aa+k,999,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]; 
  end
ldat1x=ldat1;
 for k=1:ngen
ldat1x(k,5)=0.00001;
ldat1x(k,4)=0.0000;
end
% Run prefault loadflow
if reply2 == 'n'
vt=ones(1,length(bdat1(:,1))-ngen);% prefault voltages set equal to 1
else
vt=run_prefault_loadflow(bdat1,ldat1x,econv,itermax);
end
% Short cirtuit study (Full Ybus)
%Fault at k
nl=length(ldat2(:,1));
n=length(bdat1(:,1))-ngen+length(tdat(:,1)); 
for i=1:length(bdat1(:,1))-ngen+length(tdat(:,1))
   for j=1:length(bdat1(:,1))-ngen+length(tdat(:,1))
Z2(i,j)=10^200;
	end
end
for k=1:nl
   if ldat2(k,1)>=ldat2(k,2)
      Z2(ldat2(k,2),ldat2(k,1))=complex(ldat2(k,4),ldat2(k,5));
      else
      Z2(ldat2(k,1),ldat2(k,2))=complex(ldat2(k,4),ldat2(k,5));  
end
end
for k=1:nl
   if ldat2(k,1)>=ldat2(k,2)
      su2(ldat2(k,2),ldat2(k,1))=0.5*complex(0,ldat2(k,6));
      else
      su2(ldat2(k,1),ldat2(k,2))=0.5*complex(0,ldat2(k,6));
   end
end
for i=1:length(bdat1(:,1))-ngen
   for j=1:length(bdat1(:,1))-ngen
      if i~=j
      su2(j,i)=su2(i,j);   
      end
   end
end
for i=1:length(bdat1(:,1))-ngen
   for j=1:length(bdat1(:,1))-ngen
      if i~=j
      Z2(j,i)=Z2(i,j);   
      end
   end
end
Z20=Z2;
%create Ybus
for i=1:length(bdat1(:,1))-ngen
   Yb(i,i)=0;
   for j=1:length(bdat1(:,1))-ngen
      if i~=j
      Yb(i,i)=inv(Z2(i,j))+Yb(i,i)+su2(i,j);
      Yb(i,j)=-inv(Z2(i,j));   
      end
   end
end
%Yb=Yb+sh;%Add (sh) shunt suceptance at buses
for k=1:ngen
Yb(k+length(bdat1(:,1))-2*ngen,k+length(bdat1(:,1))-2*ngen)=Yb(k+length(bdat1(:,1))-2*ngen,k+length(bdat1(:,1))-2*ngen)+inv(complex(ldat1(k,4),ldat1(k,5)));
end
Zb=inv(Yb);
Yb0=Yb;
clear U
for k=flag1+1:length(bdat1(:,1))-ngen-length(tdat(:,1));
U(:,k)=-Zb(:,k).*inv(Zb(k,k));
Ix(IM(k,2),k,k)=abs(vt(k+ngen)*Ibase*(U(IM(k,2),k)-U(k,k))/Z2(IM(k,2),k));
Ix(IM(k,1),k,k)=abs(vt(k+ngen)*Ibase*(U(IM(k,1),k)-U(k,k))/Z2(IM(k,1),k));
Ax(IM(k,2),k,k)=angle(vt(k+ngen)*Ibase*(U(IM(k,2),k)-U(k,k))/Z2(IM(k,2),k))*180/pi;
Ax(IM(k,1),k,k)=angle(vt(k+ngen)*Ibase*(U(IM(k,1),k)-U(k,k))/Z2(IM(k,1),k))*180/pi;
if IB1(k,1)==1
Ix(IB1(k,2),IB1(k,3),k)=abs(vt(k+ngen)*Ibase*((U(IB1(k,2),k)-U(IB1(k,3),k))/Z2(IB1(k,2),IB1(k,3))));  
Ax(IB1(k,2),IB1(k,3),k)=angle(vt(k+ngen)*Ibase*((U(IB1(k,2),k)-U(IB1(k,3),k))/Z2(IB1(k,2),IB1(k,3))))*180/pi;
elseif IB1(k,1)==2
Ix(IB1(k,2),IB1(k,4),k)=abs(vt(k+ngen)*Ibase*((U(IB1(k,2),k)-U(IB1(k,4),k))/Z2(IB1(k,2),IB1(k,4))));  
Ix(IB1(k,3),IB1(k,5),k)=abs(vt(k+ngen)*Ibase*((U(IB1(k,3),k)-U(IB1(k,5),k))/Z2(IB1(k,3),IB1(k,5))));
Ax(IB1(k,2),IB1(k,4),k)=angle(vt(k+ngen)*Ibase*((U(IB1(k,2),k)-U(IB1(k,4),k))/Z2(IB1(k,2),IB1(k,4))))*180/pi;  
Ax(IB1(k,3),IB1(k,5),k)=angle(vt(k+ngen)*Ibase*((U(IB1(k,3),k)-U(IB1(k,5),k))/Z2(IB1(k,3),IB1(k,5))))*180/pi;
elseif IB1(k,1)==3
Ix(IB1(k,2),IB1(k,5),k)=abs(vt(k+ngen)*Ibase*((U(IB1(k,2),k)-U(IB1(k,5),k))/Z2(IB1(k,2),IB1(k,5))));  
Ix(IB1(k,3),IB1(k,6),k)=abs(vt(k+ngen)*Ibase*((U(IB1(k,3),k)-U(IB1(k,6),k))/Z2(IB1(k,3),IB1(k,6))));  
Ix(IB1(k,4),IB1(k,7),k)=abs(vt(k+ngen)*Ibase*((U(IB1(k,4),k)-U(IB1(k,7),k))/Z2(IB1(k,4),IB1(k,7)))); 
Ax(IB1(k,2),IB1(k,5),k)=angle(vt(k+ngen)*Ibase*((U(IB1(k,2),k)-U(IB1(k,5),k))/Z2(IB1(k,2),IB1(k,5))))*180/pi;  
Ax(IB1(k,3),IB1(k,6),k)=angle(vt(k+ngen)*Ibase*((U(IB1(k,3),k)-U(IB1(k,6),k))/Z2(IB1(k,3),IB1(k,6))))*180/pi;  
Ax(IB1(k,4),IB1(k,7),k)=angle(vt(k+ngen)*Ibase*((U(IB1(k,4),k)-U(IB1(k,7),k))/Z2(IB1(k,4),IB1(k,7))))*180/pi; 
elseif IB1(k,1)==4
Ix(IB1(k,2),IB1(k,6),k)=abs(vt(k+ngen)*Ibase*((U(IB1(k,2),k)-U(IB1(k,6),k))/Z2(IB1(k,2),IB1(k,6))));  
Ix(IB1(k,3),IB1(k,7),k)=abs(vt(k+ngen)*Ibase*((U(IB1(k,3),k)-U(IB1(k,7),k))/Z2(IB1(k,3),IB1(k,7))));  
Ix(IB1(k,4),IB1(k,8),k)=abs(vt(k+ngen)*Ibase*((U(IB1(k,4),k)-U(IB1(k,8),k))/Z2(IB1(k,4),IB1(k,8)))); 
Ix(IB1(k,5),IB1(k,9),k)=abs(vt(k+ngen)*Ibase*((U(IB1(k,5),k)-U(IB1(k,9),k))/Z2(IB1(k,5),IB1(k,9)))); 
Ax(IB1(k,2),IB1(k,6),k)=angle(vt(k+ngen)*Ibase*((U(IB1(k,2),k)-U(IB1(k,6),k))/Z2(IB1(k,2),IB1(k,6))))*180/pi;  
Ax(IB1(k,3),IB1(k,7),k)=angle(vt(k+ngen)*Ibase*((U(IB1(k,3),k)-U(IB1(k,7),k))/Z2(IB1(k,3),IB1(k,7))))*180/pi;  
Ax(IB1(k,4),IB1(k,8),k)=angle(vt(k+ngen)*Ibase*((U(IB1(k,4),k)-U(IB1(k,8),k))/Z2(IB1(k,4),IB1(k,8))))*180/pi; 
Ax(IB1(k,5),IB1(k,9),k)=angle(vt(k+ngen)*Ibase*((U(IB1(k,5),k)-U(IB1(k,9),k))/Z2(IB1(k,5),IB1(k,9))))*180/pi; 
end
if IB2(k,1)==1
Ix(IB2(k,2),IB2(k,3),k)=abs(vt(k+ngen)*Ibase*((U(IB2(k,2),k)-U(IB2(k,3),k))/Z2(IB2(k,2),IB2(k,3))));  
Ax(IB2(k,2),IB2(k,3),k)=angle(vt(k+ngen)*Ibase*((U(IB2(k,2),k)-U(IB2(k,3),k))/Z2(IB2(k,2),IB2(k,3))))*180/pi;  
elseif IB2(k,1)==2
Ix(IB2(k,2),IB2(k,4),k)=abs(vt(k+ngen)*Ibase*((U(IB2(k,2),k)-U(IB2(k,4),k))/Z2(IB2(k,2),IB2(k,4))));  
Ix(IB2(k,3),IB2(k,5),k)=abs(vt(k+ngen)*Ibase*((U(IB2(k,3),k)-U(IB2(k,5),k))/Z2(IB2(k,3),IB2(k,5))));
Ax(IB2(k,2),IB2(k,4),k)=angle(vt(k+ngen)*Ibase*((U(IB2(k,2),k)-U(IB2(k,4),k))/Z2(IB2(k,2),IB2(k,4))))*180/pi;  
Ax(IB2(k,3),IB2(k,5),k)=angle(vt(k+ngen)*Ibase*((U(IB2(k,3),k)-U(IB2(k,5),k))/Z2(IB2(k,3),IB2(k,5))))*180/pi; 
elseif IB2(k,1)==3
Ix(IB2(k,2),IB2(k,5),k)=abs(vt(k+ngen)*Ibase*((U(IB2(k,2),k)-U(IB2(k,5),k))/Z2(IB2(k,2),IB2(k,5))));  
Ix(IB2(k,3),IB2(k,6),k)=abs(vt(k+ngen)*Ibase*((U(IB2(k,3),k)-U(IB2(k,6),k))/Z2(IB2(k,3),IB2(k,6))));  
Ix(IB2(k,4),IB2(k,7),k)=abs(vt(k+ngen)*Ibase*((U(IB2(k,4),k)-U(IB2(k,7),k))/Z2(IB2(k,4),IB2(k,7)))); 
Ax(IB2(k,2),IB2(k,5),k)=angle(vt(k+ngen)*Ibase*((U(IB2(k,2),k)-U(IB2(k,5),k))/Z2(IB2(k,2),IB2(k,5))))*180/pi;  
Ax(IB2(k,3),IB2(k,6),k)=angle(vt(k+ngen)*Ibase*((U(IB2(k,3),k)-U(IB2(k,6),k))/Z2(IB2(k,3),IB2(k,6))))*180/pi;  
Ax(IB2(k,4),IB2(k,7),k)=angle(vt(k+ngen)*Ibase*((U(IB2(k,4),k)-U(IB2(k,7),k))/Z2(IB2(k,4),IB2(k,7))))*180/pi;
elseif IB2(k,1)==4
Ix(IB2(k,2),IB2(k,6),k)=abs(vt(k+ngen)*Ibase*((U(IB2(k,2),k)-U(IB2(k,6),k))/Z2(IB2(k,2),IB2(k,6))));  
Ix(IB2(k,3),IB2(k,7),k)=abs(vt(k+ngen)*Ibase*((U(IB2(k,3),k)-U(IB2(k,7),k))/Z2(IB2(k,3),IB2(k,7))));  
Ix(IB2(k,4),IB2(k,8),k)=abs(vt(k+ngen)*Ibase*((U(IB2(k,4),k)-U(IB2(k,8),k))/Z2(IB2(k,4),IB2(k,8)))); 
Ix(IB2(k,5),IB2(k,9),k)=abs(vt(k+ngen)*Ibase*((U(IB2(k,5),k)-U(IB2(k,9),k))/Z2(IB2(k,5),IB2(k,9)))); 
Ax(IB2(k,2),IB2(k,6),k)=angle(vt(k+ngen)*Ibase*((U(IB2(k,2),k)-U(IB2(k,6),k))/Z2(IB2(k,2),IB2(k,6))))*180/pi;  
Ax(IB2(k,3),IB2(k,7),k)=angle(vt(k+ngen)*Ibase*((U(IB2(k,3),k)-U(IB2(k,7),k))/Z2(IB2(k,3),IB2(k,7))))*180/pi;  
Ax(IB2(k,4),IB2(k,8),k)=angle(vt(k+ngen)*Ibase*((U(IB2(k,4),k)-U(IB2(k,8),k))/Z2(IB2(k,4),IB2(k,8))))*180/pi; 
Ax(IB2(k,5),IB2(k,9),k)=angle(vt(k+ngen)*Ibase*((U(IB2(k,5),k)-U(IB2(k,9),k))/Z2(IB2(k,5),IB2(k,9))))*180/pi; 
end
end
for kk=flag1+1:length(bdat1(:,1))-ngen-length(tdat(:,1))
% Short cirtuit study (TRANSIENT CONFIG)
%Fault at k 
nl=length(ldat2(:,1));
n=length(bdat1(:,1))-ngen; 
for i=1:length(bdat1(:,1))-ngen
   for j=1:length(bdat1(:,1))-ngen
Z2(i,j)=10^200;
	end
end
for k=1:nl
   if ldat2(k,1)>=ldat2(k,2)
      Z2(ldat2(k,2),ldat2(k,1))=complex(ldat2(k,4),ldat2(k,5));
      else
      Z2(ldat2(k,1),ldat2(k,2))=complex(ldat2(k,4),ldat2(k,5));  
   end
end
  Z2(IM(kk,2),kk)=10^200;
  Z2(kk,IM(kk,2))=10^200;

for k=1:nl
   if ldat2(k,1)>=ldat2(k,2)
      su2(ldat2(k,2),ldat2(k,1))=0.5*complex(0,ldat2(k,6));
      else
      su2(ldat2(k,1),ldat2(k,2))=0.5*complex(0,ldat2(k,6));
   end
end
for i=1:length(bdat1(:,1))-ngen
   for j=1:length(bdat1(:,1))-ngen
      if i~=j
      su2(j,i)=su2(i,j);   
      end
   end
end
for i=1:length(bdat1(:,1))-ngen
   for j=1:length(bdat1(:,1))-ngen
      if i~=j
      Z2(j,i)=Z2(i,j);   
      end
   end
end
Z2x=Z2;
%create Ybus
for i=1:length(bdat1(:,1))-ngen
   Yb(i,i)=0;
   for j=1:length(bdat1(:,1))-ngen
      if i~=j
      Yb(i,i)=inv(Z2(i,j))+Yb(i,i)+su2(i,j);
      Yb(i,j)=-inv(Z2(i,j));   
      end
   end
end
%Yb=Yb+sh;%Add (sh) shunt suceptance at buses
for k=1:ngen
%Yb(k,k)=Yb(k,k)+inv(complex(ldat1(k,4),ldat1(k,5)));
Yb(k+length(bdat1(:,1))-2*ngen,k+length(bdat1(:,1))-2*ngen)=Yb(k+length(bdat1(:,1))-2*ngen,k+length(bdat1(:,1))-2*ngen)+inv(complex(ldat1(k,4),ldat1(k,5)));

end
Zb=inv(Yb);
Ybx=Yb;


   
   
U(:,kk)=-Zb(:,kk).*inv(Zb(kk,kk));
Ixx(IM(kk,1),kk,kk)=abs(Ibase*(U(IM(kk,1),kk)-U(kk,kk))/Z2(IM(kk,1),kk));
Axx(IM(kk,1),kk,kk)=angle(Ibase*(U(IM(kk,1),kk)-U(kk,kk))/Z2(IM(kk,1),kk))*180/pi;
if IB1(kk,1)==1
Ixx(IB1(kk,2),IB1(kk,3),kk)=abs(Ibase*((U(IB1(kk,2),kk)-U(IB1(kk,3),kk))/Z2(IB1(kk,2),IB1(kk,3)))); 
Axx(IB1(kk,2),IB1(kk,3),kk)=angle(Ibase*((U(IB1(kk,2),kk)-U(IB1(kk,3),kk))/Z2(IB1(kk,2),IB1(kk,3))))*180/pi;
elseif IB1(kk,1)==2
Ixx(IB1(kk,2),IB1(kk,4),kk)=abs(Ibase*((U(IB1(kk,2),kk)-U(IB1(kk,4),kk))/Z2(IB1(kk,2),IB1(kk,4))));  
Ixx(IB1(kk,3),IB1(kk,5),kk)=abs(Ibase*((U(IB1(kk,3),kk)-U(IB1(kk,5),kk))/Z2(IB1(kk,3),IB1(kk,5)))); 
Axx(IB1(kk,2),IB1(kk,4),kk)=angle(Ibase*((U(IB1(kk,2),kk)-U(IB1(kk,4),kk))/Z2(IB1(kk,2),IB1(kk,4))))*180/pi;  
Axx(IB1(kk,3),IB1(kk,5),kk)=angle(Ibase*((U(IB1(kk,3),kk)-U(IB1(kk,5),kk))/Z2(IB1(kk,3),IB1(kk,5))))*180/pi;  
elseif IB1(kk,1)==3
Ixx(IB1(kk,2),IB1(kk,5),kk)=abs(Ibase*((U(IB1(kk,2),kk)-U(IB1(kk,5),kk))/Z2(IB1(kk,2),IB1(kk,5))));  
Ixx(IB1(kk,3),IB1(kk,6),kk)=abs(Ibase*((U(IB1(kk,3),kk)-U(IB1(kk,6),kk))/Z2(IB1(kk,3),IB1(kk,6))));  
Ixx(IB1(kk,4),IB1(kk,7),kk)=abs(Ibase*((U(IB1(kk,4),kk)-U(IB1(kk,7),kk))/Z2(IB1(kk,4),IB1(kk,7)))); 
Axx(IB1(kk,2),IB1(kk,5),kk)=angle(Ibase*((U(IB1(kk,2),kk)-U(IB1(kk,5),kk))/Z2(IB1(kk,2),IB1(kk,5))))*180/pi;  
Axx(IB1(kk,3),IB1(kk,6),kk)=angle(Ibase*((U(IB1(kk,3),kk)-U(IB1(kk,6),kk))/Z2(IB1(kk,3),IB1(kk,6))))*180/pi;  
Axx(IB1(kk,4),IB1(kk,7),kk)=angle(Ibase*((U(IB1(kk,4),kk)-U(IB1(kk,7),kk))/Z2(IB1(kk,4),IB1(kk,7))))*180/pi;  
elseif IB1(kk,1)==4
Ixx(IB1(kk,2),IB1(kk,6),kk)=abs(Ibase*((U(IB1(kk,2),kk)-U(IB1(kk,6),kk))/Z2(IB1(kk,2),IB1(kk,6))));  
Ixx(IB1(kk,3),IB1(kk,7),kk)=abs(Ibase*((U(IB1(kk,3),kk)-U(IB1(kk,7),kk))/Z2(IB1(kk,3),IB1(kk,7))));  
Ixx(IB1(kk,4),IB1(kk,8),kk)=abs(Ibase*((U(IB1(kk,4),kk)-U(IB1(kk,8),kk))/Z2(IB1(kk,4),IB1(kk,8)))); 
Ixx(IB1(kk,5),IB1(kk,9),kk)=abs(Ibase*((U(IB1(kk,5),kk)-U(IB1(kk,9),kk))/Z2(IB1(kk,5),IB1(kk,9)))); 
Axx(IB1(kk,2),IB1(kk,6),kk)=angle(Ibase*((U(IB1(kk,2),kk)-U(IB1(kk,6),kk))/Z2(IB1(kk,2),IB1(kk,6))))*180/pi;  
Axx(IB1(kk,3),IB1(kk,7),kk)=angle(Ibase*((U(IB1(kk,3),kk)-U(IB1(kk,7),kk))/Z2(IB1(kk,3),IB1(kk,7))))*180/pi;  
Axx(IB1(kk,4),IB1(kk,8),kk)=angle(Ibase*((U(IB1(kk,4),kk)-U(IB1(kk,8),kk))/Z2(IB1(kk,4),IB1(kk,8))))*180/pi; 
Axx(IB1(kk,5),IB1(kk,9),kk)=angle(Ibase*((U(IB1(kk,5),kk)-U(IB1(kk,9),kk))/Z2(IB1(kk,5),IB1(kk,9))))*180/pi; 
end

% Short cirtuit study (TRANSIENT CONFIG)
%Fault at k 
nl=length(ldat2(:,1));
n=length(bdat1(:,1))-ngen; 
for i=1:length(bdat1(:,1))-ngen
   for j=1:length(bdat1(:,1))-ngen
Z2(i,j)=10^200;
	end
end
for k=1:nl
   if ldat2(k,1)>=ldat2(k,2)
      Z2(ldat2(k,2),ldat2(k,1))=complex(ldat2(k,4),ldat2(k,5));
      else
      Z2(ldat2(k,1),ldat2(k,2))=complex(ldat2(k,4),ldat2(k,5));  
   end
end
  Z2(IM(kk,1),kk)=10^200;
  Z2(kk,IM(kk,1))=10^200;

for k=1:nl
   if ldat2(k,1)>=ldat2(k,2)
      su2(ldat2(k,2),ldat2(k,1))=0.5*complex(0,ldat2(k,6));
      else
      su2(ldat2(k,1),ldat2(k,2))=0.5*complex(0,ldat2(k,6));
   end
end
for i=1:length(bdat1(:,1))-ngen
   for j=1:length(bdat1(:,1))-ngen
      if i~=j
      su2(j,i)=su2(i,j);   
      end
   end
end
for i=1:length(bdat1(:,1))-ngen
   for j=1:length(bdat1(:,1))-ngen
      if i~=j
      Z2(j,i)=Z2(i,j);   
      end
   end
end
Z2xx=Z2;
%create Ybus
for i=1:length(bdat1(:,1))-ngen
   Yb(i,i)=0;
   for j=1:length(bdat1(:,1))-ngen
      if i~=j
      Yb(i,i)=inv(Z2(i,j))+Yb(i,i)+su2(i,j);
      Yb(i,j)=-inv(Z2(i,j));   
      end
   end
end
%Y=Y+sh;%Add (sh) shunt suceptance at buses
for k=1:ngen
%Yb(k,k)=Yb(k,k)+inv(complex(ldat1(k,4),ldat1(k,5)));
Yb(k+length(bdat1(:,1))-2*ngen,k+length(bdat1(:,1))-2*ngen)=Yb(k+length(bdat1(:,1))-2*ngen,k+length(bdat1(:,1))-2*ngen)+inv(complex(ldat1(k,4),ldat1(k,5)));

end
Zb=inv(Yb);

Ybxx=Yb;
U(:,kk)=-Zb(:,kk).*inv(Zb(kk,kk));
Ixx(IM(kk,2),kk,kk)=abs(Ibase*(U(IM(kk,2),kk)-U(kk,kk))/Z2(IM(kk,2),kk)); 
Axx(IM(kk,2),kk,kk)=angle(Ibase*(U(IM(kk,2),kk)-U(kk,kk))/Z2(IM(kk,2),kk))*180/pi; 
if IB2(kk,1)==1
Ixx(IB2(kk,2),IB2(kk,3),kk)=abs(Ibase*((U(IB2(kk,2),kk)-U(IB2(kk,3),kk))/Z2(IB2(kk,2),IB2(kk,3)))); 
Axx(IB2(kk,2),IB2(kk,3),kk)=angle(Ibase*((U(IB2(kk,2),kk)-U(IB2(kk,3),kk))/Z2(IB2(kk,2),IB2(kk,3))))*180/pi; 
elseif IB2(kk,1)==2
Ixx(IB2(kk,2),IB2(kk,4),kk)=abs(Ibase*((U(IB2(kk,2),kk)-U(IB2(kk,4),kk))/Z2(IB2(kk,2),IB2(kk,4))));  
Ixx(IB2(kk,3),IB2(kk,5),kk)=abs(Ibase*((U(IB2(kk,3),kk)-U(IB2(kk,5),kk))/Z2(IB2(kk,3),IB2(kk,5))));  
Axx(IB2(kk,2),IB2(kk,4),kk)=angle(Ibase*((U(IB2(kk,2),kk)-U(IB2(kk,4),kk))/Z2(IB2(kk,2),IB2(kk,4))))*180/pi;  
Axx(IB2(kk,3),IB2(kk,5),kk)=angle(Ibase*((U(IB2(kk,3),kk)-U(IB2(kk,5),kk))/Z2(IB2(kk,3),IB2(kk,5))))*180/pi; 
elseif IB2(kk,1)==3
Ixx(IB2(kk,2),IB2(kk,5),kk)=abs(Ibase*((U(IB2(kk,2),kk)-U(IB2(kk,5),kk))/Z2(IB2(kk,2),IB2(kk,5))));  
Ixx(IB2(kk,3),IB2(kk,6),kk)=abs(Ibase*((U(IB2(kk,3),kk)-U(IB2(kk,6),kk))/Z2(IB2(kk,3),IB2(kk,6))));  
Ixx(IB2(kk,4),IB2(kk,7),kk)=abs(Ibase*((U(IB2(kk,4),kk)-U(IB2(kk,7),kk))/Z2(IB2(kk,4),IB2(kk,7))));  
Axx(IB2(kk,2),IB2(kk,5),kk)=angle(Ibase*((U(IB2(kk,2),kk)-U(IB2(kk,5),kk))/Z2(IB2(kk,2),IB2(kk,5))))*180/pi;  
Axx(IB2(kk,3),IB2(kk,6),kk)=angle(Ibase*((U(IB2(kk,3),kk)-U(IB2(kk,6),kk))/Z2(IB2(kk,3),IB2(kk,6))))*180/pi;  
Axx(IB2(kk,4),IB2(kk,7),kk)=angle(Ibase*((U(IB2(kk,4),kk)-U(IB2(kk,7),kk))/Z2(IB2(kk,4),IB2(kk,7))))*180/pi;
elseif IB2(kk,1)==4
Ixx(IB2(kk,2),IB2(kk,6),kk)=abs(Ibase*((U(IB2(kk,2),kk)-U(IB2(kk,6),kk))/Z2(IB2(kk,2),IB2(kk,6))));  
Ixx(IB2(kk,3),IB2(kk,7),kk)=abs(Ibase*((U(IB2(kk,3),kk)-U(IB2(kk,7),kk))/Z2(IB2(kk,3),IB2(kk,7))));  
Ixx(IB2(kk,4),IB2(kk,8),kk)=abs(Ibase*((U(IB2(kk,4),kk)-U(IB2(kk,8),kk))/Z2(IB2(kk,4),IB2(kk,8)))); 
Ixx(IB2(kk,5),IB2(kk,9),kk)=abs(Ibase*((U(IB2(kk,5),kk)-U(IB2(kk,9),kk))/Z2(IB2(kk,5),IB2(kk,9)))); 
Axx(IB2(kk,2),IB2(kk,6),kk)=angle(Ibase*((U(IB2(kk,2),kk)-U(IB2(kk,6),kk))/Z2(IB2(kk,2),IB2(kk,6))))*180/pi;  
Axx(IB2(kk,3),IB2(kk,7),kk)=angle(Ibase*((U(IB2(kk,3),kk)-U(IB2(kk,7),kk))/Z2(IB2(kk,3),IB2(kk,7))))*180/pi;  
Axx(IB2(kk,4),IB2(kk,8),kk)=angle(Ibase*((U(IB2(kk,4),kk)-U(IB2(kk,8),kk))/Z2(IB2(kk,4),IB2(kk,8))))*180/pi; 
Axx(IB2(kk,5),IB2(kk,9),kk)=angle(Ibase*((U(IB2(kk,5),kk)-U(IB2(kk,9),kk))/Z2(IB2(kk,5),IB2(kk,9))))*180/pi; 
end
end
cont=1;
for k=flag1+1:length(bdat1(:,1))-ngen-length(tdat(:,1))
%main nrec
for j=1:IB1(k,1)
IndexB1(cont,8)=k;
IndexB1(cont,7)=IM(k,1);
IndexB1(cont,6)=IB1(k,j+1+IB1(k,1));
IndexB1(cont,5)=IB1(k,j+1);
IndexB1(cont,4)=k-flag1;
IndexB1(cont,9)=Ix(IM(k,1),k,k);
IndexB1(cont,2)=k-flag1;
IndexB1(cont,10)=Ix(IB1(k,j+1),IB1(k,j+1+IB1(k,1)),k);
IndexB1(cont,11)=Ixx(IM(k,1),k,k);
IndexB1(cont,12)=Ixx(IB1(k,j+1),IB1(k,j+1+IB1(k,1)),k);
IndexB1(cont,17)=Ax(IM(k,1),k,k);
IndexB1(cont,18)=Ax(IB1(k,j+1),IB1(k,j+1+IB1(k,1)),k);
IndexB1(cont,19)=Axx(IM(k,1),k,k);
IndexB1(cont,20)=Axx(IB1(k,j+1),IB1(k,j+1+IB1(k,1)),k);
cont=cont+1;
end  
end
cont=1;
for k=flag1+1:length(bdat1(:,1))-ngen-length(tdat(:,1))
for j=1:IB2(k,1)
IndexB2(cont,8)=k;
IndexB2(cont,7)=IM(k,2);
IndexB2(cont,6)=IB2(k,j+1+IB2(k,1));
IndexB2(cont,5)=IB2(k,j+1);
IndexB2(cont,4)=k-flag1;
IndexB2(cont,2)=k-flag1;
IndexB2(cont,9)=Ix(IM(k,2),k,k);
IndexB2(cont,10)=Ix(IB2(k,j+1),IB2(k,j+1+IB2(k,1)),k);
IndexB2(cont,11)=Ixx(IM(k,2),k,k);
IndexB2(cont,12)=Ixx(IB2(k,j+1),IB2(k,j+1+IB2(k,1)),k);
IndexB2(cont,17)=Ax(IM(k,2),k,k);
IndexB2(cont,18)=Ax(IB2(k,j+1),IB2(k,j+1+IB2(k,1)),k);
IndexB2(cont,19)=Axx(IM(k,2),k,k);
IndexB2(cont,20)=Axx(IB2(k,j+1),IB2(k,j+1+IB2(k,1)),k);
cont=cont+1;
end
end
IndexB=vertcat(IndexB1,IndexB2);
for k=1:length(IndexB(:,1))
for j=1:length(Relay(:,2))
if Relay(j,2)==IndexB(k,5) & Relay(j,3)==IndexB(k,6)
IndexB(k,1)=Relay(j,1);
end
end
end
for k=1:length(IndexB(:,1))
for j=1:length(Relay(:,2))
if Relay(j,2)==IndexB(k,7) & Relay(j,3)==IndexB(k,8)
IndexB(k,3)=Relay(j,1);
end
end
IndexB(k,13)=k1/((IndexB(k,9)/Ip(IndexB(k,3)))^k2+k3);%Normal main
IndexB(k,14)=k1/((IndexB(k,10)/Ip(IndexB(k,1)))^k2+k3);%Normal backup
IndexB(k,15)=k1/((IndexB(k,11)/Ip(IndexB(k,3)))^k2+k3);%Trans main
IndexB(k,16)=k1/((IndexB(k,12)/Ip(IndexB(k,1)))^k2+k3);%Trans 
end
index=IndexB;  
for kw=1:length(index(:,1))
betabf(index(kw,1),index(kw,2))= index(kw,14);    
betabt(index(kw,1),index(kw,2))= index(kw,16);
betamf(index(kw,3),index(kw,4))= index(kw,13);    
betamt(index(kw,3),index(kw,4))= index(kw,15);
end 
for kw=1:length(index(:,1))/2
index(kw,22)=index(kw,3)+1;
end
for kw=length(index(:,1))/2+1:length(index(:,1))
index(kw,22)=index(kw,3)-1;
end    
for kw=1:length(index(:,1))
index(kw,23)=index(kw,4);
end
for kw=1:length(index(:,1))
index(kw,21)=(betabt(index(kw,1),index(kw,2))/betabf(index(kw,1),index(kw,2))...
-betamt(index(kw,3),index(kw,4))/betamf(index(kw,3),index(kw,4)))*betamf(index(kw,22),index(kw,4));
end
for kw=1:length(index(:,1))
index(kw,24)=-(betamt(index(kw,3),index(kw,4))/betamf(index(kw,3),index(kw,4)))*betamf(index(kw,22),index(kw,4));
end
for kw=1:length(index(:,1))
index(kw,21)=(betabt(index(kw,1),index(kw,2))/betabf(index(kw,1),index(kw,2))...
-betamt(index(kw,3),index(kw,4))/betamf(index(kw,3),index(kw,4)))*betamf(index(kw,22),index(kw,4));

index(kw,24)=(0/betabf(index(kw,1),index(kw,2))...
-betamt(index(kw,3),index(kw,4))/betamf(index(kw,3),index(kw,4)))*betamf(index(kw,22),index(kw,4));
index(kw,25)=(betabt(index(kw,1),index(kw,2))/betabf(index(kw,1),index(kw,2))...
-0/betamf(index(kw,3),index(kw,4)))*betamf(index(kw,22),index(kw,4));
end
for k=1:length(index(:,1))% null shorcircuit currents are not in reverse
if index(k,10) < 0.00001
index(k,18)=index(k,17);  
end
end 
end


