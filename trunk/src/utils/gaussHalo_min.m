function f = gaussHalo_min(p, arg1)

x  = arg1(:,1); 
y  = arg1(:,2);
dy = arg1(:,3);

A  = [ones(size(x)) gaussHaloFn(x,p(1),p(2),p(3),p(4),p(5))];
c  = A\y;
z  = A*c;
f  = norm((z-y)./dy);

              
