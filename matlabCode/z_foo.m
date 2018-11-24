function [I, firDerIn, firDerOut, secDerIn, secDerOut, secDerOutIn] = z_foo(x, y)
  uConst = 0.5;
  vConst = 0.3;
  u = x + uConst*y;
  v = x - vConst*y;
  I = u.^2 - v.^2;
  firDerIn = 2*(x + uConst*y) - 2*(x - vConst*y);
  firDerOut = 2*uConst*(x + uConst*y) + 2*vConst*(x - vConst*y);
  secDerIn = 0.0;
  secDerOut = uConst - 2*vConst*vConst;
  secDerOutIn = 0.0;
end
