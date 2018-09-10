filename = 'points-1.txt';
points = importdata(filename);
points = points';
hull = convhull(points(1,:), points(2,:));
CH = points(:, hull);
scatter(points(1, :), points(2, :), 'g')
hold on;
plot(CH(1, [1:end 1]), CH(2, [1:end 1]), 'r--o');
fprintf('convex hull area = %0.6f\n', polyarea(CH(1, :), CH(2, :)))
box = minBoundingBox(points);
plot(box(1, [1:end 1]), box(2, [1:end 1]), 'b--o')
fprintf('bounding rectangle area = %0.6f\n', polyarea(box(1, :), box(2, :)))
