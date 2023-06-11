
function diff = diffVector(triplet1, triplet2)

diff = zeros(1,3);

for position = 1:3
    if ~strcmp(triplet1(position),triplet2(position))
        diff(1,position) = 1;
    end
end

