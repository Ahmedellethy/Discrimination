% Initiate and Update The centers of the SOMs clusters

%Read the seismic event parameters data in the folowing format
% the depth of expolsions are zero
%depth	PcfHz	ScfHz	PSRcf
%08.46	10.28	08.95	1.15
%the data are tab delimited
cd 'C:\code'
formatspec= '%f%f%f%f';
% spacifiy the data file
Fld='C\Code\Training.dat'; % Change the path as approperate
%Read the data into table T
T=readtable(Fld,'Delimiter','	','Format',formatspec);
Depth=T{:,1}; 
PSRcf=T{:,2:4};
Pcf=PSRcf(:,1); %P-wave corner freq
Scf=PSRcf(:,2); %S-wave corner freq
Rcf=PSRcf(:,3); %Ratio of P-wave to S-wave corner freq
Tgt=Depth(:)~=0; %Target binary set  1 for earthquakes and 0 for explosions
Tgt=double(Tgt);
sz=size(Pcf);
TE=sz(1); %number of events
% Cluster the data using the SOM
SOMC1=SOM1(Pcf);
SOMC2=SOM2(Scf);
SOMC3=SOM3(Rcf);
SOMC4=SOM4(PSRcf);
nE1=sum(SOMC1);% number of events “hits” (nE) in each cluster
nE2=sum(SOMC2);
nE3=sum(SOMC3);
nE4=sum(SOMC4);
%Cluster support value SC=nE/TE equ. (5)
SC1=nE1/TE;
SC2=nE2/TE;
SC3=nE3/TE;
SC4=nE4/TE;
% Create arrays for Tables 3 & 4 which contains the clusters details
% No of events, no of earthquakes, no of explosions, Support, earthquakes
% confidence, and explosions confidence
C4P=zeros(6,4);
C4S=zeros(6,4);
C4R=zeros(6,4);
C9PSR=zeros(6,9);
% create arrays to hold the k-mean value of each cluster 
K4P=zeros(4,1);
K4S=zeros(4,1);
K4R=zeros(4,1);
K9=zeros(9,3);

C4P(1,:)=nE1;
C4P(4,:)=SC1;
C4S(1,:)=nE2;
C4S(4,:)=SC2;
C4R(1,:)=nE3;
C4R(4,:)=SC3;
C9PSR(1,:)=nE4;
C9PSR(4,:)=SC4;
for i=1 : 4
    Cx1=SOMC1(:,i)==1;
    K4P(i)=mean(Pcf(Cx1));% centers
    Ceq=Tgt(Cx1)==1 ;
    Cex=Tgt(Cx1)==0;
    C4P(2,i)=sum(Ceq);
    C4P(3,i)=sum(Cex);
    C4P(5,i)=C4P(2,i)/C4P(1,i); % CfEq=nEq/nE equ. (6)
    C4P(6,i)=C4P(3,i)/C4P(1,i); % CfEx=nEx/nE equ. (7)
    
    Cx1=SOMC2(:,i)==1;
    K4S(i)=mean(Scf(Cx1));% centers
    Ceq=Tgt(Cx1)==1 ;
    Cex=Tgt(Cx1)==0;
    C4S(2,i)=sum(Ceq);
    C4S(3,i)=sum(Cex);
    C4S(5,i)=C4S(2,i)/C4S(1,i); % CfEq=nEq/nE equ. (6)
    C4S(6,i)=C4S(3,i)/C4S(1,i); % CfEx=nEx/nE equ. (7)
    
    Cx1=SOMC3(:,i)==1;
    K4R(i)=mean(Rcf(Cx1));% centers
    Ceq=Tgt(Cx1)==1 ;
    Cex=Tgt(Cx1)==0;
    C4R(2,i)=sum(Ceq);
    C4R(3,i)=sum(Cex);
    C4R(5,i)=C4R(2,i)/C4R(1,i); % CfEq=nEq/nE equ. (6)
    C4R(6,i)=C4R(3,i)/C4R(1,i); % CfEx=nEx/nE equ. (7)
end 

for i=1 : 9
    Cx1=SOMC4(:,i)==1;
    K9(i,:)=mean(PSRcf(Cx1,:)); % centers of the 9 clusters
    Ceq=Tgt(Cx1)==1 ;
    Cex=Tgt(Cx1)==0;
    C9PSR(2,i)=sum(Ceq);
    C9PSR(3,i)=sum(Cex);
    C9PSR(5,i)=C9PSR(2,i)/C9PSR(1,i); % CfEq=nEq/nE equ. (6)
    C9PSR(6,i)=C9PSR(3,i)/C9PSR(1,i); % CfEx=nEx/nE equ. (7)
end
