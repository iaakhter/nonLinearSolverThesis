function p = bisection(f,gnd, Vdd, inputVolt, x0, a,b)

% provide the equation you want to solve with R.H.S = 0 form. 
% Write the L.H.S by using inline function
% Give initial guesses.
% Solves it by method of bisection.
% A very simple code. But may come handy

if f(gnd, Vdd, inputVolt, x0, a)*f(gnd, Vdd, inputVolt, x0, b)>0 
    disp('Wrong choice bro')
else
    p = (a + b)/2;
    err = abs(f(gnd, Vdd, inputVolt, x0, p));
    while err > 1e-7
      if abs(a - b) < 1e-17
        break
      end
   if f(gnd, Vdd, inputVolt, x0, a)*f(gnd, Vdd, inputVolt, x0, p)<0 
       b = p;
   else
       a = p;          
   end
    p = (a + b)/2; 
   err = abs(f(gnd, Vdd, inputVolt, x0, p));
    end
end