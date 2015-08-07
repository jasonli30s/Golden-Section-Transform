function X = ihgst2d(H,nlevel)

% Author: Jun Li. Copyright (c) 2015, All rights reserved. 
% more info@ http://goldensectiontransform.com/
% the lifting scheme of inverse 2d high golden section transform
% here matrix H is of Fn*Fn where Fn is a fibonacci number Fn>=2.

global hji; % only used by function ihgst1d(S) below.
global FBHI; % only used by function ihgst1d(S) below.

[xx,yy] = size(H);

ind = floor(log(xx*sqrt(5)+1/2)/log((sqrt(5)+1)/2)); % determine index
FBHI = filter(1,[1 -1 -1],[1 zeros(1,ind-1)]);
% FBHI = Fibonacci sequence -> [1 1 2 3 5 8...];

X=H;
for hji=nlevel:-1:1
   
   for j=1:FBHI(end-2*hji+2)
      ss = X(1:FBHI(end-2*hji),j);
      dd = X(FBHI(end-2*hji)+1:FBHI(end-2*hji+2),j);
      X(1:FBHI(end-2*hji+2),j) = ihgst1d(ss',dd')'; % inverse column transform
   end
   
   for k=1:FBHI(end-2*hji+2)
      ss = X(k,1:FBHI(end-2*hji));
      dd = X(k,FBHI(end-2*hji)+1:FBHI(end-2*hji+2));
      X(k,1:FBHI(end-2*hji+2)) = ihgst1d(ss,dd); % inverse row transform
   end
   
end


%% inverse 1d high golden section transform lifting scheme

function S = ihgst1d(ss,dd)

%% Author: Jun Li. more info@ http://goldensectiontransform.com/

global hji;
global FBHI;

N = length(ss) + length(dd);
index = 0;
g = 1;
h = 1;

hform = hword(N);

for i=1:length(hform)
   
   index = index + hform(i);
   
   if hform(i) == 2
      
      % S(index-1) = (sqrt(FBHI(2*hji-1))*ss(i)+sqrt(FBHI(2*hji))*dd(2*i-g))/sqrt(FBHI(2*hji+1));
      % S(index) = (sqrt(FBHI(2*hji))*ss(i)-sqrt(FBHI(2*hji-1))*dd(2*i-g))/sqrt(FBHI(2*hji+1));
      
      %% we can use lifting scheme here:
      
      ga = sqrt(FBHI(2*hji-1)/FBHI(2*hji));
      gb = sqrt(FBHI(2*hji-1)*FBHI(2*hji))/FBHI(2*hji+1);
      gc = sqrt(FBHI(2*hji+1)/FBHI(2*hji));
      
      d1 = gc*dd(2*i-g);
      s1 = ss(i)/gc;
      S(index) = s1 - gb*d1;
      S(index-1) = d1 + ga*S(index);
      
      g = g+1;
      
   else % hform(i) == 3
      
      % S(index-2) = (sqrt(FBHI(2*hji))*ss(i)+sqrt(FBHI(2*hji-1)/2)*dd(i+h-1))/sqrt(FBHI(2*hji+2))+dd(i+h)/sqrt(2);
      % S(index-1) = (sqrt(FBHI(2*hji-1))*ss(i)-sqrt(2*FBHI(2*hji))*dd(i+h-1))/sqrt(FBHI(2*hji+2));
      % S(index) = (sqrt(FBHI(2*hji))*ss(i)+sqrt(FBHI(2*hji-1)/2)*dd(i+h-1))/sqrt(FBHI(2*hji+2))-dd(i+h)/sqrt(2);
      
      %% we can use lifting scheme here:
      
      ha = sqrt(2*FBHI(2*hji)/FBHI(2*hji-1));
      hb = sqrt(2*FBHI(2*hji-1)*FBHI(2*hji))/FBHI(2*hji+2);
      hc = sqrt(FBHI(2*hji+2)/FBHI(2*hji-1));
      
      tdd = dd(i+h);
      td2 = hc*dd(i+h-1);
      ts2 = ss(i)/hc;
      S(index-1) = ts2 - hb*td2;
      tss = td2 + ha*S(index-1);

      td1 = sqrt(2)*tdd;
      ts1 = tss/sqrt(2);
      S(index) = ts1 - 1/2*td1;
      S(index-2) = td1 + S(index);
      
      h = h+1;
      
   end
   
end
