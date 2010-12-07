Straightening algorithm
=======================

0. Read in the image

1. Find the standard deviation of the image
2. Choose 50 pixels more than one s.d. above the mean
3. Construct complete graph where each edge is weighted wsith the distance between the points
4. Compute the Minimum Spanning Tree of this graph.
5. Find the longest path, or diameter, of the MST.
6. Construct from this an ordered sequence of control points (this will constitute the start point for the energy minimization procedure).

7. Update the control points using Eq. (10) on page 3 of the paper. Iterate until they don't move significantly.
8. Restack, using 1-pixel cutting-planes.

9. Write out the image
