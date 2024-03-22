function y = graphConnectedOrNot(x)


[m n] = size(x);
A = zeros(m,n);
for i = 1 : m
    for j = 1 : n
        if x(i, j) ~= 0
            A(i,x(i,j)) = 1;
        end
    end
end

[S C] = conncomp(A);
if all(C == ones(size(C)))
    y = 1;
    fprintf('the graph is connected!\n');
else
    y = 0;
    %fprintf('the graph is not connected!\n');
    %fprintf('result is %d\n',C);

end

end

