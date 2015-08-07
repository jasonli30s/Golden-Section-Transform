function H = lgst2d(X,nlevel)

% Author: Jun Li. Copyright (c) 2015, All rights reserved. 
% more info@ http://goldensectiontransform.com/
% the lifting scheme of 2d low golden section transform
% here matrix X is of Fn*Fn where Fn is a fibonacci number Fn>=2.

global lj; % only used by function lgst1d(S) below.
global FBL; % only used by function lgst1d(S) below.

[xx,yy] = size(X);

ind = floor(log(xx*sqrt(5)+1/2)/log((sqrt(5)+1)/2)); % determine index
FBL = filter(1,[1 -1 -1],[1 zeros(1,ind-1)]);
% FBL = Fibonacci sequence -> [1 1 2 3 5 8...];

H=X;
for lj=1:nlevel
   
   for j=1:xx
      [ss,dd] = lgst1d(H(j,1:yy)); % row transform
      H(j,1:yy) = [ss,dd];
   end
   
   for k=1:yy
      [ss,dd] = lgst1d(H(1:xx,k)'); % column transform
      H(1:xx,k) = [ss,dd]';
   end
   
   xx = FBL(end-lj); % round((sqrt(5)-1)/2*xx); 8*8 block: xx=8->5->3->2
   yy = FBL(end-lj); % round((sqrt(5)-1)/2*yy); 8*8 block: yy=8->5->3->2
   
end

%% 1d low golden section transform lifting scheme

function [ss,dd] = lgst1d(S)

%% Author: Jun Li. more info@ http://goldensectiontransform.com/

global lj;
global FBL;

index = 0;
h = 1;

lform = lword(length(S));

for i=1:length(lform)
   
   index = index + lform(i);
   
   if lform(i) == 1
      
      ss(i) = S(index);
      
   else % lform(i) == 2
      
      % ss(i) = (sqrt(FBL(lj))*S(index-1)+sqrt(FBL(lj+1))*S(index))/sqrt(FBL(lj+2));
      % dd(h) = (sqrt(FBL(lj+1))*S(index-1)-sqrt(FBL(lj))*S(index))/sqrt(FBL(lj+2));
      
      %% we can use lifting scheme here:
      
      ga = sqrt(FBL(lj)/FBL(lj+1));
      gb = sqrt(FBL(lj)*FBL(lj+1))/FBL(lj+2);
      gc = sqrt(FBL(lj+2)/FBL(lj+1));
      
      d1 = S(index-1) - ga*S(index);
      s1 = S(index) + gb*d1;
      ss(i) = gc*s1;
      dd(h) = d1/gc;
      
      h = h+1;
      
   end
   
end
