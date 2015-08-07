function hform = hword(n)

% Author: Jun Li. Copyright (c) 2015, All rights reserved. 
% more info@ http://goldensectiontransform.com/
% Type-H golden section decomposion of fibonacci number n, 
% e.g. hword(8) = [3 2 3];

if n == 2

   hform = [2];

elseif n == 3

   hform = [3];

else 

next = round((sqrt(5)-1)/2*n); 

hform = [hword(n-next),hword(next)];

end
