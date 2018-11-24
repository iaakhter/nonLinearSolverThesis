function [I, firDerIn, firDerOut, secDerIn, secDerOut, secDerOutIn] = z_foo(x, y)
  I = x.^2 - y.^2;
  firDerIn = 2*x;
  firDerOut = -2*y;
  secDerIn = 2.0;
  secDerOut = -2.0;
  secDerOutIn = 0.0;
end
