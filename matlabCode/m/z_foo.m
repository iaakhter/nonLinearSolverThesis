function z = z_foo(x, y)
  u = x + 0.5*y;
  v = x - 0.3*y;
  z = u.^2 - v.^2;
end
