function H = rhgst2d(X,nlevel)

% Author: Jun Li. Copyright (c) 2015, All rights reserved. 
% more info@ http://goldensectiontransform.com/
% the lifting scheme of 2d reverse order high golden section transform
% here matrix X is of Fn*Fn where Fn is a fibonacci number Fn>=2.
% Wythoff sequence is used here, no need to call function_hword.

global hj; % only used by function rhgst1d(S) below.
global FBH; % only used by function rhgst1d(S) below.

[xx,yy] = size(X);

ind = floor(log(xx*sqrt(5)+1/2)/log((sqrt(5)+1)/2)); % determine index
FBH = filter(1,[1 -1 -1],[1 zeros(1,ind-1)]);
% FBH = Fibonacci sequence -> [1 1 2 3 5 8...];

H=X;
for hj=1:nlevel
   
   for j=1:xx
      [ss,dd] = rhgst1d(H(j,1:yy)); % row transform
      H(j,1:yy) = [ss,dd];
   end
   
   for k=1:yy
      [ss,dd] = rhgst1d(H(1:xx,k)'); % column transform
      H(1:xx,k) = [ss,dd]';
   end
   
   xx = FBH(end-2*hj); % round((3-sqrt(5))/2*xx); 8*8 block: xx=8->3
   yy = FBH(end-2*hj); % round((3-sqrt(5))/2*yy); 8*8 block: yy=8->3
   
end

%% 1d reverse order high golden section transform lifting scheme

function [ss,dd] = rhgst1d(S)

%% Author: Jun Li. more info@ http://goldensectiontransform.com/

global hj;
global FBH;

if length(S) == 2
   
   % ss(1) = (sqrt(FBH(2*hj))*S(1)+sqrt(FBH(2*hj-1))*S(2))/sqrt(FBH(2*hj+1));
   % dd(1) = (sqrt(FBH(2*hj-1))*S(1)-sqrt(FBH(2*hj))*S(2))/sqrt(FBH(2*hj+1));
   
   %% we can use lifting scheme here:
   
   ga = sqrt(FBH(2*hj)/FBH(2*hj-1));
   gb = sqrt(FBH(2*hj)*FBH(2*hj-1))/FBH(2*hj+1);
   gc = sqrt(FBH(2*hj+1)/FBH(2*hj-1));
   
   d1 = S(1) - ga*S(2);
   s1 = S(2) + gb*d1;
   ss(1) = gc*s1;
   dd(1) = d1/gc;
   
else %length(S) == 3,5,8,13,21...
   
   for i=1:FBH(end-2*hj-1)
      
      lw = floor((1+sqrt(5))/2*i); % Lower Wythoff sequence->1,3,4,6,8...
      % uw = floor((3+sqrt(5))/2*i); % Upper Wythoff sequence->2,5,7,10,13...
      uw = lw + i;
      
      % ss(lw) = (sqrt(FBH(2*hj))*S(lw+uw-2)+sqrt(FBH(2*hj-1))*S(lw+uw-1)+sqrt(FBH(2*hj))*S(lw+uw))/sqrt(FBH(2*hj+2));
      % dd(uw-1) = (sqrt(FBH(2*hj-1))*S(lw+uw-2)-2*sqrt(FBH(2*hj))*S(lw+uw-1)+sqrt(FBH(2*hj-1))*S(lw+uw))/sqrt(2*FBH(2*hj+2));
      % dd(uw) = (S(lw+uw-2)-S(lw+uw))/sqrt(2);
      
      %% we can use lifting scheme here:
      
      td1 = S(lw+uw-2) - S(lw+uw);
      ts1 = S(lw+uw) + 1/2*td1;
      tss = sqrt(2)*ts1;
      tdd = td1/sqrt(2);
      
      ha = sqrt(2*FBH(2*hj)/FBH(2*hj-1));
      hb = sqrt(2*FBH(2*hj-1)*FBH(2*hj))/FBH(2*hj+2);
      hc = sqrt(FBH(2*hj+2)/FBH(2*hj-1));
      
      td2 = tss - ha*S(lw+uw-1);
      ts2 = S(lw+uw-1) + hb*td2;
      ss(lw) = hc*ts2;
      dd(uw-1) = td2/hc;
      dd(uw) = tdd;
      
      if (FBH(end-2*hj) ~= 1) & (i <= FBH(end-2*hj-2))
         
         % ss(uw) = (sqrt(FBH(2*hj))*S(2*uw+lw-1)+sqrt(FBH(2*hj-1))*S(2*uw+lw))/sqrt(FBH(2*hj+1));
         % dd(lw+uw) = (sqrt(FBH(2*hj-1))*S(2*uw+lw-1)-sqrt(FBH(2*hj))*S(2*uw+lw))/sqrt(FBH(2*hj+1));
         
         %% we can use lifting scheme here:
         
         ga = sqrt(FBH(2*hj)/FBH(2*hj-1));
         gb = sqrt(FBH(2*hj)*FBH(2*hj-1))/FBH(2*hj+1);
         gc = sqrt(FBH(2*hj+1)/FBH(2*hj-1));
         
         d1 = S(2*uw+lw-1) - ga*S(2*uw+lw);
         s1 = S(2*uw+lw) + gb*d1;
         ss(uw) = gc*s1;
         dd(lw+uw) = d1/gc;
         
      end
      
   end
   
end
