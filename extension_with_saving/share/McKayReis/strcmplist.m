function I = strcmplist(U,L)
% U (universe) is the set of all strings in a cell array
% L (list) is a list of the subset
% I is a vector with length of L with I(i) the position of L(i) in U.  If
% not found, then I(i) = 0

I = zeros(1,length(L));

for j = 1:length(L)
    i = find(strcmp(U,L{j}));
    if ~isempty(i)
        I(j) = i;
    end
end