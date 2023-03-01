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


def get_intervals(cell):
    # TODO this func would benefit from a check that cell has the right shape
    """

    Undo what get_corners does to a list of intervals,
    e.g., recover the intervals [ [a,b], [c,d], [e,f] ]
    from the corners [ [b,d,f], [a,c,e] ].
    Remarkabbly, the operation is involutive.
    
    """
    foo = np.array(cell)[:,::-1]
    return foo.T



def measure_cell( cell, volume=True ):
    ##### ~~~
    ## ~~~~ WARNING: as written, this function can return negative length and volume
    #### ~~~
    if not (np.shape(cell)[0]==2 and len(np.shape(cell))==2):
        raise ValueError("Cell not formatted correctly.")
    upper_right  =  np.array( cell[0], dtype=np.float64 )
    lower_left   =  np.array( cell[1], dtype=np.float64 )
    lengths  =  (upper_right - lower_left).tolist()
    if volume:
        Lebesgue_measure  =  np.prod( lengths, dtype=np.float64 )
    else:
        Lebesgue_measure  =  None
    return lengths, Lebesgue_measure


def measure_multiple_cells( list_of_cells, volume=True ):
        lengths_of_each_side  =  list()
        Lebesgue_measure  =  list()
        for cell in list_of_cells:
            lengths, volume  =  measure_cell(cell,volume)
            lengths_of_each_side.append(lengths)
            Lebesgue_measure.append(volume)
        return (
            np.array( lengths_of_each_side, dtype=np.float64 ),
            Lebesgue_measure
            )


class ntangle_beta:
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
    #   define the fields
    #
    def __init__(self,list_of_lists_of_intervals):
        #
        #~~~ first, convert each list of intervals into a pair of corners
        list_of_cells = [
                get_corners(list_of_intervals).tolist()
                for list_of_intervals in list_of_lists_of_intervals
            ]
        #
        # ~~~ check a couple of simple size requirements
        shape = np.shape(list_of_cells)
        if not len(shape)==3:
            raise ValueError("Cells nor formatted correctly.")
        if not shape[1]==2:
            raise ValueError("Cells nor formatted correctly.")
        #
        # ~~~ gather data
        lengths_of_each_side, Lebesgue_measure  =  measure_multiple_cells(list_of_cells)
        weights  =  Lebesgue_measure / np.sum(Lebesgue_measure)
        list_of_cells = np.array( list_of_cells, dtype=np.float64 )
        #
        # ~~~ throw away and degenerate cells
        non_degenerate_ones = np.where( np.all(lengths_of_each_side>0, axis=0) )
        if len(non_degenerate_ones[0])<len(list_of_cells):
            list_of_cells = list_of_cells[non_degenerate_ones]
            warnings.warn("Degenerate cells discarded.")
        #
        # ~~~ assign attributes
        self.cells  =  list_of_cells
        self.intervals  =  np.array( list_of_lists_of_intervals, dtype=np.float64 )
        self.number_of_cells  =  len(list_of_cells)
        self.lengths_of_each_side  =  lengths_of_each_side
        self.Lebesgue_measure  =  Lebesgue_measure
        self.weights  =  weights
        self.dimension  =  len(list_of_cells[0][0])
        #
        # ~~~ run a sanity check or two
        self.check_nondegeneracy()
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
                    raise ValueError("Cells have non-disjoint interiors, i.e., they overlap. This will cause issues when sampling.")
    #
    #   sample a point at random from a randomly chosen cell
    #
    def sample( self, N, with_time=None ):
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
        if with_time==None:
            pass # do nothing
        else:
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
Omega = ntangle_beta(good)

# (ex.5) provokes ValueError: Cells have non-disjoint interiors...
Prompts_Error = ntangle_beta(bad)

# (ex.6) equivalent to ex.1, returns (False, None)
open_cells_intersect( Omega.cells[0], Omega.cells[1] )

# (ex.7) sample 10 points chosen uniformlt at random from Omega = ([0,5]x[0,2]) \setminus ([0,1)X[0,1))
Omega.sample(10)


#
