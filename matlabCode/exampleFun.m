function [I, firDerIn, firDerOut, secDerIn, secDerOut] = exampleFun(in, out)
  I = (in - out)*(in - out);
  firDerIn = 2*(in - out);
  firDerOut = -2*(in - out);
  secDerIn = 2.0;
  secDerOut = 2.0;
end
