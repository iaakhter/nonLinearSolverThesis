function handles = saddle(z)
  if(nargin<1)
    z = @(x,y)(z_default(x,y));
  end
  %close(1);
  x = -0.8:0.1:1.2;
  y = -0.3:0.1:1.0;
  xx = repmat(x', 1, length(y));
  yy = repmat(y, length(x), 1);
  zz = z(xx, yy);
  figure(1);
  surf(xx, yy, zz);
  hold('on');
  x1 = zeros(4, length(x)-1);
  y1 = zeros(4, length(x)-1);
  for i = 1:size(x1,2)
    x1(:,i) = [x(i); x(i+1); x(i+1); x(i)];
    y1(:,i) = [y(1); y(1);   y(end); y(end)];
  end
  z1 = z(x1, y1);
  c1 = mean(z1, 1);
  p1 = patch(x1, y1, z1, c1);
  x2 = zeros(4, length(y)-1);
  y2 = zeros(4, length(y)-1);
  for i = 1:size(y2,2)
    x2(:,i) = [x(1); x(1);   x(end); x(end)];
    y2(:,i) = [y(i); y(i+1); y(i+1); y(i)];
  end
  z2 = z(x2, y2);
  c2 = mean(z2, 1);
  p2 = patch(x2, y2, z2, c2);
  alpha(p1, 0.3);
  alpha(p2, 0.3);
  handles = {p1, p2};
  c3 = zeros(4,1,3);
  for i=1:size(c3,1)
    for j = 1:size(c3,2)
      for k = 1:size(c3,3)
        c3(i,j,k) = my_color(i,j,k);
      end
    end
  end
  p3 = patch(x2(:,4), y2(:,4), z2(:,4), c3);
end

function z = z_default(x, y)
  z = x.^2 - y.^2;
end

function c = my_color(i,j,k)
  if(k==1) c = 1.0;
  else c = 0.0;
  end
end
