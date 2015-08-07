function H = rlgst2d(X,nlevel)

% Author: Jun Li. Copyright (c) 2015, All rights reserved. 
% more info@ http://goldensectiontransform.com/
% the lifting scheme of 2d reverse order low golden section transform
% here matrix X is of Fn*Fn where Fn is a fibonacci number Fn>=2.
% Wythoff sequence is used here, no need to call function_lword.

global lj; % only used by function rlgst1d(S) below.
global FBL; % only used by function rlgst1d(S) below.

[xx,yy] = size(X);

ind = floor(log(xx*sqrt(5)+1/2)/log((sqrt(5)+1)/2)); % determine index
FBL = filter(1,[1 -1 -1],[1 zeros(1,ind-1)]);
% FBL = Fibonacci sequence -> [1 1 2 3 5 8...];

H=X;
for lj=1:nlevel
   
   for j=1:xx
      [ss,dd] = rlgst1d(H(j,1:yy)); % row transform
      H(j,1:yy) = [ss,dd];
   end
   
   for k=1:yy
      [ss,dd] = rlgst1d(H(1:xx,k)'); % column transform
      H(1:xx,k) = [ss,dd]';
   end
   
   xx = FBL(end-lj); % round((sqrt(5)-1)/2*xx); 8*8 block: xx=8->5->3->2
   yy = FBL(end-lj); % round((sqrt(5)-1)/2*yy); 8*8 block: yy=8->5->3->2
   
end

%% 1d reverse order low golden section transform lifting scheme

function [ss,dd] = rlgst1d(S)

%% Author: Jun Li. more info@ http://goldensectiontransform.com/

global lj;
global FBL;

for i=1:FBL(end-lj-1)
   
   lw = floor((1+sqrt(5))/2*i); % Lower Wythoff sequence->1,3,4,6,8...
   % uw = floor((3+sqrt(5))/2*i); % Upper Wythoff sequence->2,5,7,10,13...
   uw = lw + i;
   
   % ss(lw) = (sqrt(FBL(lj+1))*S(uw-1)+sqrt(FBL(lj))*S(uw))/sqrt(FBL(lj+2));
   % dd(i) = (sqrt(FBL(lj))*S(uw-1)-sqrt(FBL(lj+1))*S(uw))/sqrt(FBL(lj+2));
   
   %% we can use lifting scheme here:
   
   ga = sqrt(FBL(lj+1)/FBL(lj));
   gb = sqrt(FBL(lj)*FBL(lj+1))/FBL(lj+2);
   gc = sqrt(FBL(lj+2)/FBL(lj));
   
   d1 = S(uw-1) - ga*S(uw);
   s1 = S(uw) + gb*d1;
   ss(lw) = gc*s1;
   dd(i) = d1/gc;
   
   if (FBL(end-lj) ~= 1) & (i <= FBL(end-lj-2))
      
      ss(uw) = S(lw+uw);
      
   end
   
end
