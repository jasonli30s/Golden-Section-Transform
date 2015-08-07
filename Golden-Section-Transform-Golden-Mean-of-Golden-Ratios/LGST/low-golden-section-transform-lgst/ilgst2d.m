function X = ilgst2d(H,nlevel)

% Author: Jun Li. Copyright (c) 2015, All rights reserved. 
% more info@ http://goldensectiontransform.com/
% the lifting scheme of inverse 2d low golden section transform
% here matrix H is of Fn*Fn where Fn is a fibonacci number Fn>=2.

global lji; % only used by function ilgst1d(S) below.
global FBLI; % only used by function ilgst1d(S) below.

[xx,yy] = size(H);

ind = floor(log(xx*sqrt(5)+1/2)/log((sqrt(5)+1)/2)); % determine index
FBLI = filter(1,[1 -1 -1],[1 zeros(1,ind-1)]);
% FBLI = Fibonacci sequence -> [1 1 2 3 5 8...];

X=H;
for lji=nlevel:-1:1
   
   for j=1:FBLI(end-lji+1)
      ss = X(1:FBLI(end-lji),j);
      dd = X(FBLI(end-lji)+1:FBLI(end-lji+1),j);
      X(1:FBLI(end-lji+1),j) = ilgst1d(ss',dd')'; % inverse column transform
   end
   
   for k=1:FBLI(end-lji+1)
      ss = X(k,1:FBLI(end-lji));
      dd = X(k,FBLI(end-lji)+1:FBLI(end-lji+1));
      X(k,1:FBLI(end-lji+1)) = ilgst1d(ss,dd); % inverse row transform
   end
   
end


%% inverse 1d low golden section transform lifting scheme

function S = ilgst1d(ss,dd)

%% Author: Jun Li. more info@ http://goldensectiontransform.com/

global lji;
global FBLI;

N = length(ss) + length(dd);
index = 0; 
h = 1;

lform = lword(N); 

for i=1:length(lform)
   
   index = index + lform(i);
   
   if lform(i) == 1
      
      S(index) = ss(i);
      
   else % lform(i) == 2
      
      % S(index-1) = (sqrt(FBLI(lji))*ss(i)+sqrt(FBLI(lji+1))*dd(h))/sqrt(FBLI(lji+2));
      % S(index) = (sqrt(FBLI(lji+1))*ss(i)-sqrt(FBLI(lji))*dd(h))/sqrt(FBLI(lji+2));
      
      %% we can use lifting scheme here:
      
      ga = sqrt(FBLI(lji)/FBLI(lji+1));
      gb = sqrt(FBLI(lji)*FBLI(lji+1))/FBLI(lji+2);
      gc = sqrt(FBLI(lji+2)/FBLI(lji+1));
      
      d1 = gc*dd(h);
      s1 = ss(i)/gc;
      S(index) = s1 - gb*d1;
      S(index-1) = d1 + ga*S(index);
      
      h = h+1;
      
   end
   
end
