function [yX] = CMnet(CMin)
%CMnet
%
% [yX] = CMnet(CMin) takes these arguments:
%  CMin = Qx12 matrix
%   x = Qx4 matrix, input #1
% and returns:
%   yX = Qx3 matrix, output #1
% where Q is the number of samples.

x1=double(CMin(:,1:4));

% ===== NEURAL NETWORK CONSTANTS =====

% Input 1
x1_step1.xoffset = [0;0;0;0];
x1_step1.gain = [2;2;2;2];
x1_step1.ymin = -1;

% Layer 1
b1 = [2.3906175785346790441;-1.9835659269036667318;-1.3522823159076751143;-0.85395894188781817924;-0.5101854374381671331;0.34553619179808275907;-0.84612598216325618505;1.3721819584682415361;-1.9551486808735329959;2.6913440525311642659];
IW1_1 = [-1.1183976398223216542 -1.7779324186694938081 0.64760710499288598463 1.3623714239614508337;1.6800798390627782464 -0.22864642003981996976 -0.63996434623183129364 -1.6890471957200414721;0.14680540474352546387 2.4585863852939544572 -0.17592753191615151387 0.49092921030252739678;1.5590797667854565667 -1.1771778056665467282 0.045179409615837454184 -1.5282205261690902098;1.7022399694693375327 0.50797861375666408801 -0.5927466877100688869 -1.5926231281091030034;1.655269586820511396 -1.6965804405527260545 0.10949692161726767847 -0.70102362344438684527;-1.3320022461651082057 -2.0171567208808340865 0.0038908476050679235113 -0.56574303785384483234;0.58052288229685178322 -1.2499035198029293525 1.9352283183109546538 -0.7317818945596428204;-0.76676889295913730482 -2.0621873981822580468 -1.1064953104865098421 -0.17039250369079317338;1.7090138089525093168 -0.20449215932967867992 -1.0995195004039701736 1.0478906795887283909];

% Layer 2
b2 = -0.58768454987404861178;
LW2_1 = [0.30390235391903369644 0.33992947588812200133 0.63240915560251709238 -0.30403087458978067525 -0.44044909531680415116 0.18821449973601253602 -0.14705524503286565574 0.69881743444431076995 -0.14903574762162100087 -0.33237857353853178521];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 2;
y1_step1.xoffset = 0;

% ===== SIMULATION ========

% Dimensions
Q = size(x1,1); % samples

% Input 1
x1 = x1';
xp1 = mapminmax_apply(x1,x1_step1);

% Layer 1
a1 = tansig_apply(repmat(b1,1,Q) + IW1_1*xp1);

% Layer 2
a2 = repmat(b2,1,Q) + LW2_1*a1;

% Output 1
y1 = mapminmax_reverse(a2,y1_step1);
y1 = int8(y1');
sc=CMin(:,[5,7,9,11]);
z1=CMin(:,[6,8,10,12]);
yX=zeros(Q,3);
yX(:,1)=y1;
 for i=1:Q
    for k=1 : 4
         if CMin(i,k) == y1(i)
            yX(i,2)=max([sc(i,k),yX(i,2)]);
            yX(i,3)=max([z1(i,k),yX(i,3)]);
         end
    end
 end
 
end

% ===== MODULE FUNCTIONS ========

% Map Minimum and Maximum Input Processing Function
function y = mapminmax_apply(x,settings)
y = bsxfun(@minus,x,settings.xoffset);
y = bsxfun(@times,y,settings.gain);
y = bsxfun(@plus,y,settings.ymin);
end

% Sigmoid Symmetric Transfer Function
function a = tansig_apply(n,~)
a = 2 ./ (1 + exp(-2*n)) - 1;
end

% Map Minimum and Maximum Output Reverse-Processing Function
function x = mapminmax_reverse(y,settings)
x = bsxfun(@minus,y,settings.ymin);
x = bsxfun(@rdivide,x,settings.gain);
x = bsxfun(@plus,x,settings.xoffset);
end

