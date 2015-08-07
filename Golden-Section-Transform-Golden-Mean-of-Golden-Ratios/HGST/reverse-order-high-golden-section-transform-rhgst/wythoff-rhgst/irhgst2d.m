function X = irhgst2d(H,nlevel)

% Author: Jun Li. Copyright (c) 2015, All rights reserved. 
% more info@ http://goldensectiontransform.com/
% the lifting scheme of inverse 2d reverse order high golden section transform
% here matrix H is of Fn*Fn where Fn is a fibonacci number Fn>=2.
% Wythoff sequence is used here, no need to call function_hword.

global hji; % only used by function irhgst1d(S) below.
global FBHI; % only used by function irhgst1d(S) below.

[xx,yy] = size(H);

ind = floor(log(xx*sqrt(5)+1/2)/log((sqrt(5)+1)/2)); % determine index
FBHI = filter(1,[1 -1 -1],[1 zeros(1,ind-1)]);
% FBHI = Fibonacci sequence -> [1 1 2 3 5 8...];

X=H;
for hji=nlevel:-1:1
   
   for j=1:FBHI(end-2*hji+2)
      ss = X(1:FBHI(end-2*hji),j);
      dd = X(FBHI(end-2*hji)+1:FBHI(end-2*hji+2),j);
      X(1:FBHI(end-2*hji+2),j) = irhgst1d(ss',dd')'; % inverse column transform
   end
   
   for k=1:FBHI(end-2*hji+2)
      ss = X(k,1:FBHI(end-2*hji));
      dd = X(k,FBHI(end-2*hji)+1:FBHI(end-2*hji+2));
      X(k,1:FBHI(end-2*hji+2)) = irhgst1d(ss,dd); % inverse row transform
   end
   
end


%% inverse 1d reverse order high golden section transform lifting scheme

function S = irhgst1d(ss,dd)

%% Author: Jun Li. more info@ http://goldensectiontransform.com/

global hji;
global FBHI;

if (length(ss)+length(dd)) == 2
   
   % S(1) = (sqrt(FBHI(2*hji))*ss(1)+sqrt(FBHI(2*hji-1))*dd(1))/sqrt(FBHI(2*hji+1));
   % S(2) = (sqrt(FBHI(2*hji-1))*ss(1)-sqrt(FBHI(2*hji))*dd(1))/sqrt(FBHI(2*hji+1));
   
   %% we can use lifting scheme here:
   
   ga = sqrt(FBHI(2*hji)/FBHI(2*hji-1));
   gb = sqrt(FBHI(2*hji)*FBHI(2*hji-1))/FBHI(2*hji+1);
   gc = sqrt(FBHI(2*hji+1)/FBHI(2*hji-1));
   
   d1 = gc*dd(1);
   s1 = ss(1)/gc;
   S(2) = s1 - gb*d1;
   S(1) = d1 + ga*S(2);


else % (length(ss)+length(dd)) == 3,5,8,13,21...
   
   for i=1:FBHI(end-2*hji-1)
      
      lw = floor((1+sqrt(5))/2*i); % Lower Wythoff sequence->1,3,4,6,8...
      % uw = floor((3+sqrt(5))/2*i); % Upper Wythoff sequence->2,5,7,10,13...
      uw = lw + i;
      
      % S(lw+uw-2) = sqrt(FBHI(2*hji))/sqrt(FBHI(2*hji+2))*ss(lw)+sqrt(FBHI(2*hji-1))/sqrt(2*FBHI(2*hji+2))*dd(uw-1)+dd(uw)/sqrt(2);
      % S(lw+uw-1) = (sqrt(FBHI(2*hji-1))*ss(lw)-sqrt(2*FBHI(2*hji))*dd(uw-1))/sqrt(FBHI(2*hji+2));
      % S(lw+uw) = sqrt(FBHI(2*hji))/sqrt(FBHI(2*hji+2))*ss(lw)+sqrt(FBHI(2*hji-1))/sqrt(2*FBHI(2*hji+2))*dd(uw-1)-dd(uw)/sqrt(2);
      
      %% we can use lifting scheme here:
      
      ha = sqrt(2*FBHI(2*hji)/FBHI(2*hji-1));
      hb = sqrt(2*FBHI(2*hji-1)*FBHI(2*hji))/FBHI(2*hji+2);
      hc = sqrt(FBHI(2*hji+2)/FBHI(2*hji-1));
      
      tdd = dd(uw);
      td2 = hc*dd(uw-1);
      ts2 = ss(lw)/hc;
      S(lw+uw-1) = ts2 - hb*td2;
      tss = td2 + ha*S(lw+uw-1);
      
      td1 = sqrt(2)*tdd;
      ts1 = tss/sqrt(2);
      S(lw+uw) = ts1 - 1/2*td1;
      S(lw+uw-2) = td1 + S(lw+uw);
      
      if (FBHI(end-2*hji) ~= 1) & (i <= FBHI(end-2*hji-2))
         
         % S(2*uw+lw-1) = (sqrt(FBHI(2*hji))*ss(uw)+sqrt(FBHI(2*hji-1))*dd(lw+uw))/sqrt(FBHI(2*hji+1));
         % S(2*uw+lw) = (sqrt(FBHI(2*hji-1))*ss(uw)-sqrt(FBHI(2*hji))*dd(lw+uw))/sqrt(FBHI(2*hji+1));
         
         %% we can use lifting scheme here:
         
         ga = sqrt(FBHI(2*hji)/FBHI(2*hji-1));
         gb = sqrt(FBHI(2*hji)*FBHI(2*hji-1))/FBHI(2*hji+1);
         gc = sqrt(FBHI(2*hji+1)/FBHI(2*hji-1));
         
         d1 = gc*dd(lw+uw);
         s1 = ss(uw)/gc;
         S(2*uw+lw) = s1 - gb*d1;
         S(2*uw+lw-1) = d1 + ga*S(2*uw+lw);
         
      end
      
   end
   
end
