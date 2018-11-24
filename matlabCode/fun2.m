%author: Mark Greenstreet

function fun(a, b, u)
  if(nargin < 1) a = 5; end
  if(nargin < 2) b = 0.0; end
  if(nargin < 3) u = 4; end
  if(numel(u) == 1) u = [-u, u]; end
  % I_y = tanh(-a*x)+b-y;  (I_y == 0) => y = tanh(-a*x) + b
  % I_x = tanh(-a*y)-b-x;  (I_x == 0) => x = tanh(-a*y) - b
  uu = u(1):((u(2)-u(1))/5000.0):u(2);
  p1 = plot(uu, tanh(-a*uu) + b, 'b-');  % I_y == 0
  hold('on');
  p2 = plot(tanh(-a*uu) - b, uu, 'r-');  % I_x == 0
  hold('off');
  legend([p1, p2], {'I_y == 0', 'I_x == 0'});
end
