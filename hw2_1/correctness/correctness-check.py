#!/usr/bin/env python3
from itertools import groupby
import numpy as np
import argparse
import os

def get_frames ( out_file ):
    """
    Read the out file produces from hw2's -o option and convert it to
    a list of frames. Each frame is a list of the particles in the
    frame. Each particle is described by an (x, y) pair.
    """

    out_frames = []

    # Read Header
    num_parts, box_size = next(out_file).split()
    num_parts, box_size = int(num_parts), float(box_size)

    # Read Body
    file_sections = groupby(out_file, lambda x: x and not x.isspace())
    frame_sections = (x[1] for x in file_sections if x[0])
    for frame_section in frame_sections:
        out_frames.append( [] )
        for line in frame_section:
            x, y = line.split()
            x, y = float(x), float(y)
            out_frames[-1].append( (x, y) )
        assert( len( out_frames[-1] ) == num_parts )

    return out_frames


def calculate_dist ( frames1, frames2 ):
    """
    Calculate a list of average euclidean distances between
    corresponding points in corresponding frames. This function returns
    a list where element i is the average distances between particles
    in frame i of frames1 and frame i of frames2.
    """

    assert( len( frames1 ) == len( frames2 ) )

    avg_dists = []

    for i, frames in enumerate( zip( frames1, frames2 ) ):
        frame1, frame2 = frames
        test_xs, test_ys = list( zip( *frame1 ) )
        verf_xs, verf_ys = list( zip( *frame2 ) )
        test_xs, test_ys = np.array( test_xs ), np.array( test_ys )
        verf_xs, verf_ys = np.array( verf_xs ), np.array( verf_ys )
        dists = np.sqrt(   np.square( test_xs - verf_xs )
                         + np.square( test_ys - verf_ys ) )
        avg_dists.append( np.mean( dists ) )

    return avg_dists


def check_conditions ( avg_dists ):
    """
    Check correctness conditions for your code's output.

    1) The average distances in the first 50 steps must be less than 1e-7.

    2) The average distances over the whole 1000 steps must be less than 0.025.
    """

    assert( np.mean( avg_dists[:50] ) < 3e-7 )
    assert( np.mean( avg_dists ) < 0.025 )


if __name__ == "__main__":
    description = "Check the correctness of HW2 output."
    parser = argparse.ArgumentParser( description = description )
    parser.add_argument( "test_out", type = argparse.FileType(),
                         help = "Your code's output" )
    parser.add_argument( "verf_out", type = argparse.FileType(),
                         help = "O(n^2) code's output" )
    args = parser.parse_args()

    test_frames = get_frames( args.test_out )
    verf_frames = get_frames( args.verf_out )
    avg_dists = calculate_dist( test_frames, verf_frames )
    check_conditions( avg_dists )
    test_out_name = os.path.basename(args.test_out.name)
    verf_out_name = os.path.basename(args.verf_out.name)

    print(f'Check between {test_out_name} and {verf_out_name} succeeded!')
