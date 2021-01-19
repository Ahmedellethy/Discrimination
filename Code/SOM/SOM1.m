function [y1] = SOM1(x1)
%MYNEURALNETWORKFUNCTION neural network simulation function.
%
% Auto-generated by MATLAB, 15-Jan-2021 21:14:36.
% 
% [y1] = myNeuralNetworkFunction(x1) takes these arguments:
%   x = Qx1 matrix, input #1
% and returns:
%   y = Qx4 matrix, output #1
% where Q is the number of samples.

%#ok<*RPMT0>

% ===== NEURAL NETWORK CONSTANTS =====

% Layer 1
IW1_1 = [4.1540517241379300017;7.5576870748299329605;5.8081437125748500705;10.501456310679612827];

% ===== SIMULATION ========

% Input 1
x1 = x1';
% no processing

% Layer 1
z1 = negdist_apply(IW1_1,x1);
a1 = compet_apply(z1);

% Output 1
y1 = a1;
y1 = y1';
end

% ===== MODULE FUNCTIONS ========

% Negative Distance Weight Function
function z = negdist_apply(w,p,~)
  [S,R] = size(w);
  Q = size(p,2);
  if isa(w,'gpuArray')
    z = iNegDistApplyGPU(w,p,R,S,Q);
  else
    z = iNegDistApplyCPU(w,p,S,Q);
  end
end
function z = iNegDistApplyCPU(w,p,S,Q)
  z = zeros(S,Q);
  if (Q<S)
    pt = p';
    for q=1:Q
      z(:,q) = sum(bsxfun(@minus,w,pt(q,:)).^2,2);
    end
  else
    wt = w';
    for i=1:S
      z(i,:) = sum(bsxfun(@minus,wt(:,i),p).^2,1);
    end
  end
  z = -sqrt(z);
end
function z = iNegDistApplyGPU(w,p,R,S,Q)
  p = reshape(p,1,R,Q);
  sd = arrayfun(@iNegDistApplyGPUHelper,w,p);
  z = -sqrt(reshape(sum(sd,2),S,Q));
end
function sd = iNegDistApplyGPUHelper(w,p)
  sd = (w-p) .^ 2;
end

% Competitive Transfer Function
function a = compet_apply(n,~)
  if isempty(n)
    a = n;
  else
    [S,Q] = size(n);
    nanInd = any(isnan(n),1);
    a = zeros(S,Q,'like',n);
    [~,maxRows] = max(n,[],1);
    onesInd = maxRows + S*(0:(Q-1));
    a(onesInd) = 1;
    a(:,nanInd) = NaN;
  end
end
