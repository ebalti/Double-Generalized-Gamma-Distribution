function D = Delta(a,b)
v = zeros(1,a);
for ii=1:a
v(ii) = (b+ii-1)/a;    
end
D = v;
end