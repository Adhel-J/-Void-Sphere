# Void-Sphere
# Modifications made
* Added void Sphere object
* Changes made to InterceptRaySphere to check if intersections are inside a void sphere and skip if they are
* Changes made to InterceptRayTriangle to skip intersection that are inside a void sphere
* Changes made to ClosestIntersection to skip rendering of void objects and move intersections that are withing a void sphere to move them to just outside the void sphere
* Added function to check if object is inside a void sphere
* Added function to check if intersection is inside a void sphere