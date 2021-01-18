% Machine Discriminator (MD) script
% Written by: Ahmed Lethy
% alethy@nriag.sci.eg         ahmedellethy@gmail.com
% National Research Institute of Astronomy and Geophysics (NRIAG), Helwan, Egypt.
%Read the seismic event parameters data in the folowing format
%PcfHz	ScfHz	PSRcf
%10.28	08.95	1.15
%the data are tab delimited
cd 'C:\Code'
formatspec= '%f%f%f';
% spacifiy the data file
Fld='C\Code\test-data.dat'; % Change the path as approperate
%Read the data into table T
Tt=readtable(Fld,'Delimiter','	','Format',formatspec); 
PSRcf=Tt{:,:};
Pcf=PSRcf(:,1); %P-wave corner freq
Scf=PSRcf(:,2); %S-wave corner freq
Rcf=PSRcf(:,3); %Ratio of P-wave to S-wave corner freq
sz=size(Pcf);
TE=sz(1); %number of events
%This data is feed to the ANNs presented in Figure 4 and listed in Table 2.
ANN1=net1(Pcf')';
ANN2=net2(Pcf')';
ANN3=net3(Scf')';
ANN4=net4(Scf')';
ANN5=net5(Rcf')';
ANN6=net6(Rcf')';
ANN7=net7(PSRcf')';
ANN8=net8(PSRcf')';
% the results are combined using Eq. (4). 
ANNC1=ANN1>2 & ANN2>0.5;
ANNC2=ANN3>2 & ANN4>0.5;
ANNC3=ANN5>2 & ANN6>0.5;
ANNC4=ANN7>2 & ANN8>0.5;
% finding the holding clusters
CMin=zeros(TE,12);
CMin(:,1)=ANNC1;
CMin(:,2)=ANNC2;
CMin(:,3)=ANNC3;
CMin(:,4)=ANNC4;

for i= 1 : TE
        dd=pdist2(K4P,Pcf(i,:));
        [xm,ik]=min(dd);
        CMin(i,5)=C4P(4,ik);% support measure
        if ANNC1(i)==1
            CMin(i,6)=C4P(5,ik); % confidence of earthquake
        else
            CMin(i,6)=C4P(6,ik); % confidence of explosion
        end
    
        dd=pdist2(K4S,Scf(i,:));
        [xm,ik]=min(dd);
        CMin(i,7)=C4S(4,ik);% support measure
        if ANNC2(i)==1
            CMin(i,8)=C4S(5,ik); % confidence of earthquake
        else
            CMin(i,8)=C4S(6,ik); % confidence of explosion
        end
    
        dd=pdist2(K4R,Rcf(i,:));
        [xm,ik]=min(dd);
        CMin(i,9)=C4R(4,ik);% support measure
        if ANNC3(i)==1
            CMin(i,10)=C4R(5,ik); % confidence of earthquake
        else
            CMin(i,10)=C4R(6,ik); % confidence of explosion
        end
        
        dd=pdist2(K9,PSRcf(i,:));
        [xm,ik]=min(dd);
        CMin(i,11)=C9PSR(4,ik);% support measure
        if ANNC4(i)==1
            CMin(i,12)=C9PSR(5,ik); % confidence of earthquake
        else
            CMin(i,12)=C9PSR(6,ik); % confidence of explosion
        end
end
CMout=CMnet(CMin);
CMout(:,2:3)=CMout(:,2:3)*100; % present the support and confidence as pecentage
% Save the results to CM.txt file
% The results are the event type in the 1st col where 1 for earthquakes and 0 for explosions
% The Support value in the 2nd col and the confidence are in the 3rd col
fid = fopen('CM.txt','w');
    for row=1:TE
        fprintf(fid,'%1u\t%4.2f\t%4.2f\r\n',CMout(row,:));
    end
fclose(fid);