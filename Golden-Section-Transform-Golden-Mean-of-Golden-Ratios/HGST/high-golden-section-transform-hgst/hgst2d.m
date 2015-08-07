function H = hgst2d(X,nlevel)

% Author: Jun Li. Copyright (c) 2015, All rights reserved. 
% more info@ http://goldensectiontransform.com/
% the lifting scheme of 2d high golden section transform
% here matrix X is of Fn*Fn where Fn is a fibonacci number Fn>=2.

global hj; % only used by function hgst1d(S) below.
global FBH; % only used by function hgst1d(S) below.

[xx,yy] = size(X);

ind = floor(log(xx*sqrt(5)+1/2)/log((sqrt(5)+1)/2)); % determine index
FBH = filter(1,[1 -1 -1],[1 zeros(1,ind-1)]);
% FBH = Fibonacci sequence -> [1 1 2 3 5 8...];

H=X;
for hj=1:nlevel
   
   for j=1:xx
      [ss,dd] = hgst1d(H(j,1:yy)); % row transform
      H(j,1:yy) = [ss,dd];
   end
   
   for k=1:yy
      [ss,dd] = hgst1d(H(1:xx,k)'); % column transform
      H(1:xx,k) = [ss,dd]';
   end
   
   xx = FBH(end-2*hj); % round((3-sqrt(5))/2*xx); 8*8 block: xx=8->3
   yy = FBH(end-2*hj); % round((3-sqrt(5))/2*yy); 8*8 block: yy=8->3
   
end

%% 1d high golden section transform lifting scheme

function [ss,dd] = hgst1d(S)

%% Author: Jun Li. more info@ http://goldensectiontransform.com/

global hj;
global FBH;

index = 0;
g = 1;
h = 1;

hform = hword(length(S));

for i=1:length(hform)
   
   index = index + hform(i);
   
   if hform(i) == 2
      
      % ss(i) = (sqrt(FBH(2*hj-1))*S(index-1)+sqrt(FBH(2*hj))*S(index))/sqrt(FBH(2*hj+1));
      % dd(2*i-g) = (sqrt(FBH(2*hj))*S(index-1)-sqrt(FBH(2*hj-1))*S(index))/sqrt(FBH(2*hj+1));
      
      %% we can use lifting scheme here:
      
      ga = sqrt(FBH(2*hj-1)/FBH(2*hj));
      gb = sqrt(FBH(2*hj-1)*FBH(2*hj))/FBH(2*hj+1);
      gc = sqrt(FBH(2*hj+1)/FBH(2*hj));
      
      d1 = S(index-1) - ga*S(index);
      s1 = S(index) + gb*d1;
      ss(i) = gc*s1;
      dd(2*i-g) = d1/gc;
      
      g = g+1;
      
   else % hform(i) == 3
      
      % ss(i) = (sqrt(FBH(2*hj))*S(index-2)+sqrt(FBH(2*hj-1))*S(index-1)+sqrt(FBH(2*hj))*S(index))/sqrt(FBH(2*hj+2));
      % dd(i+h-1) = (sqrt(FBH(2*hj-1))*S(index-2)-2*sqrt(FBH(2*hj))*S(index-1)+sqrt(FBH(2*hj-1))*S(index))/sqrt(2*FBH(2*hj+2));
      % dd(i+h) = (S(index-2)-S(index))/sqrt(2);
      
      %% we can use lifting scheme here:
      
      td1 = S(index-2) - S(index);
      ts1 = S(index) + 1/2*td1;
      tss = sqrt(2)*ts1;
      tdd = td1/sqrt(2);
      
      ha = sqrt(2*FBH(2*hj)/FBH(2*hj-1));
      hb = sqrt(2*FBH(2*hj-1)*FBH(2*hj))/FBH(2*hj+2);
      hc = sqrt(FBH(2*hj+2)/FBH(2*hj-1));
            
      td2 = tss - ha*S(index-1);
      ts2 = S(index-1) + hb*td2;
      ss(i) = hc*ts2;
      dd(i+h-1) = td2/hc;
      dd(i+h) = tdd;
      
      h = h+1;
      
   end
   
end
