% ODOCR - Dataset pm.deoliveira@uniandes.edu.co
% ---------------------------------------
% Braga, A. S., and J. Tome Saraiva. 
% "Coordination of overcurrent directional relays in meshed networks using the Simplex method." 
% Proceedings of 8th Mediterranean Electrotechnical Conference on Industrial Applications in Power Systems, 
% Computer Science and Telecommunications (MELECON 96). Vol. 3. IEEE, 1996.
% 8-bus test system 7-lines 2-generator/transformers
%
%            G1
%            |
%          ----- N14
%            |
%            T1
%            |           
%       N1 -----  N2 -----  N3 ---- 
%           |||  L1   | |  L2  | |
%          12||1__N7_2| |3_N8_4| 5
%           | \___N13__          |
%          N12         \14       N9
%        L6 |  __N11__ | _N10__  |L3
%          11 10  L5  9|8  L4  7 6
%           | |       |||      | |
%       N6 -----  N5 ----- N4 ----- 
%                     |
%                     T2
%                     |
%                   ----- N15
%                     |
%                     G2
% lines K k=1,2,3,4,5,6,7
% nodes B b=1,...,14
% relays R i=1,2,3,4,5,6,7,8,9,10,11,12,13,14
%        Backup	Primary	Backup	Primary		
%        Ri	Lk	Ri	Lk	Bs	Br	Bs  Br	
% Index=[11	1	1	1	6	12	1	7	;
%        14	1	1	1	5	13	1	7	;
%        1	2	3	2	1	7	2	8	;
%        3	3	5	3	2	8	3	9	;
%        5	4	7	4	3	9	4	10	;
%        7	5	9	5	4	10	5	11	;
%        13	5	9	5	1	13	5	11	;
%        9	6	11	6	5	11	6	12	;
%        2	7	13	7	2	7	1	13	;
%        11	7	13	7	6	12	1	13	;
%        4	1	2	1	3	8	2	7	;
%        6	2	4	2	4	9	3	8	;
%        8	3	6	3	5	10	4	9	;
%        10	4	8	4	6	11	5	10	;
%        13	4	8	4	1	13	5	10	;
%        12	5	10	5	1	12	6	11	;
%        2	6	12	6	2	7	1	12	;
%        14	6	12	6	5	13	1	12	;
%        7	7	14	7	4	10	5	13	;
%        10	7	14	7	6	11	5	13	];
global Dmin Co  nr bdat ldat k1 k2 k3 Ip Sbase Vbase Zbase Ibase neval econv itermax tdat nlf  
%% Relay data %%%%%%
% Linear solution (No transient network configurations)
% Solution method::
% Urdaneta, Alberto J., Ramon Nadira, and LG Perez Jimenez. 
% "Optimal coordination of directional overcurrent relays in interconnected power systems." 
% IEEE Transactions on Power Delivery 3.3 (1988): 903-911.
% Solution also reported in:
% Sorrentino, Elmer, and José Vicente Rodríguez.
%"A novel and simpler way to include transient configurations in optimal coordination of directional overcurrent protections."
% Electric Power Systems Research 180 (2020): 106127.
trad=1;
Co=.3;% allowable coordination interval (seconds)
% %Relay curve settings, Standard Inverse (SI) IEC 60255
 %k1=.14;k2=0.02;k3=-1;
%Relay curve settings, Very Inverse (VI) IEC 60255
k1=13.5;k2=1;k3=-1;
%Pickup currents (kA)    
Ip=[0.5;0.12;0.24;0.12;0.12;0.24;0.12;0.6;0.18;0.12	;0.12;0.18;0.12;0.12];% 
Dmin=0.1; % minimal dial setting
nr=length(Ip); 
dictionodes=[1 2 3 4 5 6 7 8   9 10 11 12 13 14 15; 
             1 3 4 5 6 2 9 10 11 12 13 14 15  7 8];   % original case numbering     
dictiorelays=[1 2 3 4 5 6 7 8 9 10 11 12 13 14; 
              2 9 3 10 4 11 5 12 6 13 1 8 14 7];  % original case numbering 
% Very inverse
 D=[0.383000000000000;0.157000000000000;0.536000000000000;0.613000000000000;0.612000000000000;0.541000000000000;0.157000000000000;0.313000000000000;0.634000000000000;0.351000000000000;0.351000000000000;0.635000000000000;0.724000000000000;0.723000000000000];
%% System Data
Sbase=150;%MVA
Vbase=150;%kV
Vbaselow=10;%kV
Zbase=Vbase^2/Sbase;%ohms
Zbaselow=Vbaselow^2/Sbase;%ohms
Ibase=Sbase/(sqrt(3)*Vbase);%kA
econv=0.00000001;
itermax=100;
% Generator rated powers and reactances    
Rt=[0 0];% pu
Xt=[.04 .04];% pu
Xg=[.15 .15];%*Zbaselow/Zbase;% pu
Rg=[0 0];%pu
ngen=length(Xg);
% lines
Le=[100 70 80 100 110 90 100]; 
R=[.004 .0057 .005 .005 .0045 .0044 .005].*Le/Zbase; %line resistances (ohms)
X=[ .05 .0714 .0563 .045 .0409 .05 .05].*Le/Zbase; %line reactance (ohms)
B=[0 0 0 0 0 0 0 0]*Zbase; %line total susceptance (siemens)
nlf=length(R);
ldat=[  1 3 10001 Rg(1) Xg(1) 0   0   0;
        2 7 10002 Rg(2) Xg(2) 0   0   0;
        3 4 90002 R(2) X(2) B(2)  0   0;
        4 5 90003 R(3) X(3) B(3)  0   0;
        5 6 90004 R(4) X(4) B(4)  0   0;        
        6 7 90005 R(5) X(5) B(5)  0   0;
        7 8 90005 R(6) X(6) B(6)  0   0; 
        8 3 90001 R(1) X(1) B(1)  0   0;
        3 7 90006 R(7) X(7) B(7)  0   0;  ];
tdat=[ max(ldat(:,1))+1 ldat(1,2) 90031 Rt(1) Xt(1) 0  0   0;
       max(ldat(:,1))+2 ldat(2,2) 90036 Rt(2) Xt(2) 0  0   0];
bdat=[  1 1001  1 0    0 0.0000 0.0000  1   0    0   0  0   0     0  0    1.0  1.0   0  0      0 0  0 0;
        2 1001  2 150/150 0 0.0000 0.0000 1 0    0   0  0   0     0  0     .9  1.10  0  0      0 0  0 0;
        3 9001  3 0    0 0.0000 0.0000  1   0    0   0  0   0     0  0     .9  1.10  0  0      0 0  0 0;
        4 9002  3 0    0 60/150 40/150  1   0    0   0  0   0     0  0     .9  1.10  0  0      0 0  0 0;
        5 9003  3 0    0 70/150 40/150  1   0    0   0  0   0     0  0     .9  1.10  0  0      0 0  0 0;
        6 9001  3 0    0 70/150 50/150  1   0    0   0  0   0     0  0     .9  1.10  0  0      0 0  0 0;
        7 9002  3 0    0 0.0000 0.0000  1   0    0   0  0   0     0  0     .9  1.10  0  0      0 0  0 0;
        8 9003  3 0    0 40/150 20/150  1   0    0   0  0   0     0  0     .9  1.10  0  0      0 0  0 0;];     
