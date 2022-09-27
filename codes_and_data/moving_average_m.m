function y = moving_average_m(x, M)
N = size(x,1);

h = (1/M) * ones(M, 1);
y = zeros(size(x));

for i = 1:N
   avg =  conv(x(i,:), h, 'same');
   y(i,ceil(M/2) : size(x,2)+1-ceil(M/2)) = avg(ceil(M/2) : size(x,2)+1-ceil(M/2));
end