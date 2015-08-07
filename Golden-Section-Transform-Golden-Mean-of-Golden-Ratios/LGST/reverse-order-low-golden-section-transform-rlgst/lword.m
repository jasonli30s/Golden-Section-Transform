function lform = lword(n)

% Author: Jun Li. Copyright (c) 2015, All rights reserved. 
% more info@ http://goldensectiontransform.com/
% Type-L golden section decomposion of fibonacci number n, 
% e.g. lword(8) = [1 2 2 1 2];

if n == 1

   lform = [1];

elseif n == 2

   lform = [2];

else 

next = round((sqrt(5)-1)/2*n); 

lform = [lword(n-next),lword(next)];

end
