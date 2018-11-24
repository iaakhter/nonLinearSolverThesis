% author: Mark Greenstreet

function fun(a, b, u) 
  if(nargin < 1) a = 0.3; end
  if(nargin < 2) b = 0.1; end
  if(nargin < 3) u = 4; end
  if(numel(u) == 1) u = [-u, u]; end
  x = u(1):((u(2)-u(1))/4000.0):u(2);
  plot(x, tanh(x), 'b-');
  hold('on');
  plot(u, a*u+b, 'r-');
end
