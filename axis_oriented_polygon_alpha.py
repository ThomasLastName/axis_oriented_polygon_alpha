import numpy as np
import warnings


def get_corners(intervals):
    # TODO this func would benefit from a check that intervals has the right shape
    """
    
    A cell in R^d can be described, either, by its d interval sides, or by 2 d-dimensional corners.
    For example, the cell [a,b]x[c,d]x[e,f] in R^3 could be described by the two corners [ (a,c,e), (b,d,f) ].
    This function converts a list of intervals like [ [a,b], [c,d], [e,f] ] into corners [ [b,d,f], [a,c,e] ].
    The transopose alone gives [ [a,c,e], [b,d,f] ].
    Therefore, to get the "upper right" corner on top (an arbitrary convention), we first reverse each row.
   
    """
    foo = np.array(intervals)[:,::-1]
    return foo.T


def get_intervals(corners):
    # TODO this func would benefit from a check that cell has the right shape
    """
    Undo what get_corners does to a list of intervals,
    e.g., recover the intervals [ [a,b], [c,d], [e,f] ]
    from the corners [ [b,d,f], [a,c,e] ].
    Remarkabbly, the operation is involutive.
    
    """
    foo = np.array(corners)[:,::-1]
    return foo.T



def measure_cell( corners, volume=True ):
    # todo: is it possible to include a check that data is supplied in corner format, or possible to accept data supplied in either format?
    # todo: would it be good to add a required "dimension" argument in order to force the user to be careful?
    ##### ~~~
    ## ~~~~ WARNING: as written, this function can return negative length and volume
    #### ~~~
    if not (np.shape(corners)[0]==2 and len(np.shape(corners))==2):
        raise ValueError("Cell not formatted correctly.")
    upper_right  =  np.array( corners[0], dtype=np.float64 )
    lower_left   =  np.array( corners[1], dtype=np.float64 )
    lengths  =  (upper_right - lower_left).tolist()
    if volume:
        Lebesgue_measure  =  np.prod( lengths, dtype=np.float64 )
    else:
        Lebesgue_measure  =  None
    return lengths, Lebesgue_measure


def measure_multiple_cells( list_of_corner_pairs, volume=True ):
    lengths_of_each_side  =  list()
    Lebesgue_measures  =  list()
    for corners in list_of_corner_pairs:
        lengths, volume  =  measure_cell(corners,volume)
        lengths_of_each_side.append(lengths)
        Lebesgue_measures.append(volume)
    return (
        np.array( lengths_of_each_side, dtype=np.float64 ),
        np.array( Lebesgue_measures, dtype=np.float64 )
        )


class axis_oriented_polygon:
    '''
    fields:
        cells,
        intervalss,
        number_of_cells,
        lengths_of_each_side,
        Lebesgue_measure,
        weights,
        dimension
    '''
    #
    #   define the attributes
    #
    def __init__( self, list_of_lists_of_intervals, enforce_disjoint_interiors=True ):
        # todo: add a required argument format (one of "intervals" or "corners") with no default value thus forcing the user to be careful
        # todo: in order to further force the user to be careful, don't infer "dimension" but rather force the user to specify it as a third argument
        #
        #~~~ first, convert each list of intervals into a pair of corners
        list_of_cells = np.array([
                get_corners(list_of_intervals).tolist()
                for list_of_intervals in list_of_lists_of_intervals
            ], dtype=np.float64 )
        list_of_lists_of_intervals  =  np.array( list_of_lists_of_intervals, dtype=np.float64 )
        #
        # ~~~ check a couple of simple size requirements
        shape = np.shape(list_of_cells)
        if not len(shape)==3:
            raise ValueError("Cells not formatted correctly.")
        if not shape[1]==2:
            raise ValueError("Cells not formatted correctly.")
        #
        # ~~~ gather data
        lengths_of_each_side, Lebesgue_measure  =  measure_multiple_cells(list_of_cells)
        weights  =  Lebesgue_measure / np.sum(Lebesgue_measure)
        #
        # ~~~ now, we have enough data to detect if there are degenerate cells and discard them
        non_degenerate_ones  =  np.where( np.all(lengths_of_each_side>0, axis=1) )
        if len(non_degenerate_ones[0])<len(list_of_cells):
            warnings.warn("Degenerate cells detected and discarded.")
            list_of_cells  =  list_of_cells[non_degenerate_ones]
            list_of_lists_of_intervals  =  list_of_lists_of_intervals[non_degenerate_ones]
            lengths_of_each_side  =  lengths_of_each_side[non_degenerate_ones]
            Lebesgue_measure  =  Lebesgue_measure[non_degenerate_ones]
            weights  =  weights[non_degenerate_ones]
        #
        # ~~~ assign attributes
        self.cells  =  list_of_cells
        self.intervals  =  list_of_lists_of_intervals
        self.number_of_cells  =  len(list_of_cells)
        self.lengths_of_each_side  =  lengths_of_each_side
        self.Lebesgue_measure  =  Lebesgue_measure
        self.weights  =  weights
        self.dimension  =  len(list_of_cells[0][0])
        #
        # ~~~ run a sanity check or two
        # self.check_nondegeneracy()
        if enforce_disjoint_interiors:
            self.check_disjointness()
    #
    #   Check that each specified cell is non-empty
    #
    def check_nondegeneracy(self):
        non_degenerate_ones = np.where( np.all(self.lengths_of_each_side>0, axis=0) )
        if len(non_degenerate_ones[0])<len(self.cells):
            warnings.warn("Degenerate cells detected.")
    #
    #   Check that the specified cells do not overlap
    #
    def check_disjointness(self):
        # todo rewrite to get_intervals
        # computational burden is O(n^2)
        n = self.number_of_cells
        for i in range(n):
            for j in range(i):
                # see comments within the function open_cells_intersect
                cell = np.array(self.cells[i]).T
                CELL = np.array(self.cells[j]).T
                intervals = cell[:,::-1]    # reverse each row
                INTERVALS = CELL[:,::-1]    # reverse each row
                a = intervals[:, 0]     # vecotor of lower bounds of
                b = intervals[:, 1]     # vector of upper bounds
                A = INTERVALS[:, 0]
                B = INTERVALS[:, 1]
                if np.all( (B > a) & (b > A) ):
                    warnings.warn("Cells have non-disjoint interiors, i.e., they overlap. This will cause issues when sampling.")
    #
    #   sample a point at random from a randomly chosen cell
    #
    def sample( self, N, with_time=None ):
        # todo: rewrite to allow with_time to be a numeric value, too
        #
        # ~~~ decide which cell to sample
        number_of_cells  =  len(self.cells)
        which_cell  =  np.random.choice(np.arange(number_of_cells), N, p=self.weights)
        #
        # ~~~ sample the appropriate number of points from each cell
        lengths  =  np.array(self.lengths_of_each_side)
        low_corners  =  np.array([cell[1] for cell in self.cells])
        sampled_points  =  np.random.uniform( size=(N,self.dimension) )
        sampled_points  =  sampled_points * lengths[which_cell] + low_corners[which_cell]
        #
        #~~~ add a time component if so desired (same as in the parent class)
        if with_time is None:
            pass # do nothing
        else:
            if isinstance( with_time, int ):
                with_time  =  [with_time, with_time]
            if len(with_time)==2:
                initial_time  =  with_time[0]
                final_time    =  with_time[1]
                times = np.random.uniform( initial_time, final_time, size=(N,1) )
                sampled_points = np.hstack(( sampled_points, times ))
        return sampled_points


def open_cells_intersect( cell, CELL ):
    # todo rewrite to accept cells of either possible format, and call get_intervals instead of repeating its guts
    """
    first, convert the list of corners into a list of intervals,
    e.g., convert corners=[ [b,d,f], [a,c,e] ] into intervals [ [a,b], [c,d], [e,f] ]
    note that the transopose alone gives [ [b,a], [d,c], [f,e] ]
    """
    cell = np.array(cell).T
    CELL = np.array(CELL).T
    intervals = cell[:,::-1]    # reverse each row
    INTERVALS = CELL[:,::-1]    # reverse each row
    # Extract lower and upper bounds for both interval arrays
    a = intervals[:, 0]     # vecotor of lower bounds of
    b = intervals[:, 1]     # vector of upper bounds
    A = INTERVALS[:, 0]
    B = INTERVALS[:, 1]
    # return a,b,A,B        # for debugging
    if False:               # why is this not equivalent???
        a = np.array(cell[1])
        b = np.array(cell[0])
        A = np.array(CELL[1])
        B = np.array(CELL[0])
    """
    Claim. Assuming a<A and b<B, we have that...
       ... intervals (a,b) and (A,B) intersect IFF both B>a and b>A.
    
    Pf.
        Indeed, in general (a,b)\cap(A,B) = (max{a,A},min{b,B})
            (see https://math.stackexchange.com/q/545035)
        where (x,y) is interpreted as empty IFF y \leq x.
        Thus, the intersection is non-empty if and only if max{a,A}< min{b,B},
        in other words if and only if both max{a,A}<b and max{a,A}<B.
        However, a<b by assumption, so the former holds IFF A<b,
        As well, A<B by assumption, so the latter holds IFF a<B.
        By substitution, (a,b)\cap(A,B) is non-empty IFF both A<b and a<B.
            QED.
    This claim tells us how to compute whether open intervals I_j and H_j intersect.
    Think of cell=\prod_jI_j and CELL=\prod_jH_j.
    Using the claim, we now compute the indices j for which I_j \cap H_j is non-empty.
    """
    indices_with_intersection = (B > a) & (b > A)
    """
    Finally, recall the distributive law
        (\prod_jI_j) \cap (\prod_jH_j) = \prod_j (I_j \cap H_j)
    and that a Cartesian product is non-empty IFF so is each "Cartesian multiplicand."
    Therfore, cell \cap CELL is non-empty IFF I_j \cap H_j is non-empty for every j,
    that is to say IFF every index j is in indices_with_intersection.
    """
    non_empty = np.all(indices_with_intersection)   # True IFF indices_with_intersection is all True
    if non_empty:
        intersection = np.vstack(( np.maximum(A,a), np.minimum(B,b) )).tolist()
        intersection = get_intervals(intersection)
    else:
        intersection = None
    return non_empty, intersection




class boundary(axis_oriented_polygon):
    def __init__( self, list_of_lists_of_intervals, list_of_integers, list_of_floats ):
        super().__init__( list_of_lists_of_intervals=list_of_lists_of_intervals, enforce_disjoint_interiors=False )
        self.missing_dimension  =  np.array(list_of_integers)
        self.missing_value  =  np.array(list_of_floats)
    def sample( self, N, with_time=None ):
        #
        # ~~~ decide which cell to sample (same as in the parent class)
        d  =  self.dimension
        number_of_cells  =  len(self.cells)
        which_cell  =  np.random.choice(np.arange(number_of_cells), N, p=self.weights)
        #
        # ~~~ sample the appropriate number of points from each cell (same as in the parent class)
        lengths  =  self.lengths_of_each_side
        low_corners  =  np.array([cell[1] for cell in self.cells])
        sampled_points  =  np.random.uniform( size=(N,d) )
        sampled_points  =  sampled_points * lengths[which_cell] + low_corners[which_cell]
        #
        #~~~ embed into one higher spatial dimension
        """
        Transform every row of the array sampled_points such that:
            - the i-th row of sampled_points is replaced by a vector of length d+1
            - where the new entry in this row is derived by taking the integer value k=which_cell[i]
            - and the new vector of length d+1 has value missing_value[k] in entry missing_dimension[k].
        For example:

# given the data...
d = 4
x_row = np.zeros(d)
k = 2
val = 1
x_new = np.insert(x_row, k, val)
x_new

# ...the command...
np.insert( A.flatten(), np.arange(10)*4 + missing_dimension, missing_value ).reshape((10,5))

# ...creates the same matrix as the one called "A_new" in the block:
A = np.ones((10,4))
A_new = np.zeros((10,5))
missing_dimension = np.arange(10)%5
missing_value = np.arange(10)+10
for row in range(A.shape[0]):
    k = missing_dimension[row]
    val = missing_value[row]
    A_new[row,:] = np.insert(A[row,:], k, val)

        """
        sampled_points  =  np.insert(
                sampled_points.flatten(),
                np.arange(N)*d + self.missing_dimension[which_cell],
                self.missing_value[which_cell]
            ).reshape((N,d+1))
        #
        #~~~ add a time component if so desired (same as in the parent class)
        if with_time is None:
            pass # do nothing
        else:
            if isinstance( with_time, int ):
                with_time  =  [with_time, with_time]
            if len(with_time)==2:
                initial_time  =  with_time[0]
                final_time    =  with_time[1]
                times = np.random.uniform( initial_time, final_time, size=(N,1) )
                sampled_points = np.hstack(( sampled_points, times ))
        return sampled_points



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
            [ [1,5], [0,2] ]    # the same as above
    ]


# (ex.1) returns (False, None)
open_cells_intersect( get_corners(good[0]), get_corners(good[1]) )


# (ex.2) returns
#(True, array([[1, 2],
#       [1, 2]]))
# because the cells overlap on [1,2]x[1,2]
open_cells_intersect( get_corners(bad[0]), get_corners(bad[1]) )


# (ex.3) returns (False, None)
open_cells_intersect( get_corners(bigger[0]), get_corners(bigger[1]) )


# (ex.4) runs and doen't return anything
Omega = axis_oriented_polygon(good)


# (ex.5) provokes the UserWarning: Cells have non-disjoint interiors...
prompts_warning = axis_oriented_polygon(bad)


# (ex.6) provokes the UserWarning: Degenerate cells detected and discarded.
prompts_warning = axis_oriented_polygon(very_bad)
# returns
#array([[[1., 5.],
#        [0., 2.]]])
prompts_warning.intervals

# (ex.7) equivalent to ex.1, returns (False, None)
open_cells_intersect( Omega.cells[0], Omega.cells[1] )


# (ex.8) sample 10 points chosen uniformly at random from Omega = ([0,5]x[0,2]) \setminus ([0,1)X[0,1))
Omega.sample(10)


# (ex.9) boundary of the cell [-9,0]x[1,2]
bd = boundary( [
            [[1,2]],
            [[1,2]],
            [[-9,0]],
            [[-9,0]]
        ],
        [ 0,   0,   1,   1],
        [ 0., -9.,  1.,  2.]
    )
bd.sample( 20, with_time=[10,11] )


