function myshuffledMatrix = sub_shuffleAdjacencyMatrix(matrix)
%This function takes an adjacency matrix and shuffles the weights between
%the nodes to produce a random connectivity matrix
% Adham Elshahabi 3rd july 2013

myshuffledMatrix = zeros(size(matrix));
carrieriVector = [];
rng('shuffle')

for i = 1:size(matrix,2)
    carrieriVector = [carrieriVector matrix(i,[1:i-1])];
end
shuffledVector = carrieriVector(randperm(length(carrieriVector)));

for j=2:size(matrix,2)
    
    myshuffledMatrix(j,[1:j-1]) = shuffledVector(sum(1:j-2)+1:sum(1:j-2)+j-1);
    myshuffledMatrix([1:j-1],j) = shuffledVector(sum(1:j-2)+1:sum(1:j-2)+j-1);
end


end