function [d Q mgQ] = sampson_dist(pt, parms )

% pt = point [x y]
% parms = ellipse parameters
x = pt(:,1);
y = pt(:,2);

[a b c d e f] = deal(parms(1),parms(2),parms(3),parms(4),parms(5),parms(6));

Q = a*x.^2 + b*x.*y + c*y.^2 + d*x + e*y + f;

mgQ = sqrt( (4*a^2 + b^2)*x.^2 + (b^2 + 4*c^2)*y.^2 + 4*b*( a+c )*x.*y + ...
    (4*a*d + 2*b*e).*x + (4*c*e + 2*b*d)*y + d^2 + e^2 );
d = Q./mgQ;