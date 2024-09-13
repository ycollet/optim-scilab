// 3-bar structure

// Nodal coordinate
// [x1, y1, z1; x2, y2, z2; x3, y3, z3;....]
//
//  1    2     3
// +-----+-----+


nodes=[-1,   0, 0; ... // Node 1
       0,    0, 0; ... // Node 2
        1,   0, 0];    // Node 3

// Cross section A and Material Property (Young's modulus E)
A = .001; // m^2
E = 2E11;  // Pa

// node 1, node 2,area, E
members = [1,2,A,E; ...
           2,3,A,E];


// Non-active nodes are fixed.
active_nodes = [2];
loads = [2,0,0,-6000];

// dof = 1    for a heat problem
// dof = 2, 3 for a truss problem (2D or 3D)
// dof = 6    for a frame problem

dof = 3;

// Calculate nodal displacement vector
// Also graphic display will be invoked.

displacements = solve_fem_displ(nodes,members,active_nodes,loads,dof);

plot_fem(nodes,members,displacements,dof);
member_forces = forces_truss(nodes,members,displacements);

