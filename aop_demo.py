import axis_oriented_polygon_alpha as aop


good = [
            [ [0,1], [1,2] ],   # the Cartesian product [0,1]x[1,2]
            [ [1,5], [0,2] ]    # the Cartesian product [1,5]x[0,2]
    ]
bigger = [
            [ [-1,1], [-1,3] ], # the Cartesian product [-1,1]x[-1,3]
            [ [1,5], [0,2] ]    # the Cartesian product [1,5]x[0,2]
    ]
bad = [
            [ [0,2], [1,2] ],   # the Cartesian product [0,2]x[1,2]
            [ [1,5], [0,2] ]    # the Cartesian product [1,5]x[0,2]
    ]
very_bad = [
            [ [2,0], [1,2] ],   # non-sense
            [ [2,0], [1,2] ],   # non-sense
            [ [1,5], [0,2] ]    # the same as above
    ]


# (ex.1) returns (False, None)
aop.open_cells_intersect( aop.get_corners(good[0]), aop.get_corners(good[1]) )


# (ex.2) returns
#(True, array([[1, 2],
#       [1, 2]]))
# because the cells overlap on [1,2]x[1,2]
aop.open_cells_intersect( aop.get_corners(bad[0]), aop.get_corners(bad[1]) )


# (ex.3) returns (False, None)
aop.open_cells_intersect( aop.get_corners(bigger[0]), aop.get_corners(bigger[1]) )


# (ex.4) runs and doen't return anything
Omega = aop.axis_oriented_polygon(good)


# (ex.5) provokes the UserWarning: Cells have non-disjoint interiors...
prompts_warning = aop.axis_oriented_polygon(bad)


# (ex.6) provokes the UserWarning: Degenerate cells detected and discarded.
prompts_warning = aop.axis_oriented_polygon(very_bad)
# returns
#array([[[1., 5.],
#        [0., 2.]]])
prompts_warning.intervals

# (ex.7) equivalent to ex.1, returns (False, None)
aop.open_cells_intersect( Omega.cells[0], Omega.cells[1] )


# (ex.8) sample 10 points chosen uniformly at random from Omega = ([0,5]x[0,2]) \setminus ([0,1)X[0,1))
Omega.sample(10)


# (ex.9) boundary of the cell [-9,0]x[1,2]
bd = aop.boundary( [
            [[1,2]],
            [[1,2]],
            [[-9,0]],
            [[-9,0]]
        ],
        [ 0,   0,   1,   1],
        [ 0., -9.,  1.,  2.]
    )
bd.sample( 20, with_time=[10,11] )
