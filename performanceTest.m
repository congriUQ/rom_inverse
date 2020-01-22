%performance tests
M = rand(16, 10000);
l1 = diag(1:16);
l2 = (1:16)';
N = 10000;

tic
for i = 1:N
    tmp = diag(l2)*M;
end
t1 = toc

tic
for i = 1:N
    tmp = l2.*M;
end
t2 = toc