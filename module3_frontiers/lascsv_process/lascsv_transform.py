#!/usr/bin/env python3
import optparse
import numpy
import pandas
import shapefile
import os
from shapely.geometry.polygon import Polygon

def rotate_mat(xy, centroid, theta, center = True):
    """
    Rotate points given theta.

    Parameters
    ----------
    xy : numpy.ndarray
        A position matrix of shape (2,p).
        Where:
            The first row is X positions.
            The second row is Y positions.
            'p' is the number of points.
    centroid : numpy.ndarray
        A center matrix of shape (2,1) or (2,p). Points are rotated around
        this (these) points.
    theta : float
        Angle to rotate points by in radians. **NOT degrees**

    Returns
    -------
    xy_prime : numpy.ndarray
        A rotated matrix of shape (2,p).
        Where:
            The first row is transformed X positions.
            The second row is transformed Y positions.
            'p' is the number of points.
    """
    # get sin(theta), cos(theta)
    sin_theta = numpy.sin(theta)
    cos_theta = numpy.cos(theta)

    # construct rotation matrix
    R = numpy.array([[cos_theta, -sin_theta],
                     [sin_theta,  cos_theta]])

    # rotation of points
    # Algebraic representation:
    # x_new = x*cos(theta) - y*sin(theta)
    # y_new = x*sin(theta) + y*cos(theta)
    # Matrix representation:
    # [[x_new], = [[cos(theta), -sin(theta)], x [[x],
    #  [y_new]]    [sin(theta),  cos(theta)]]    [y]]
    xy_prime = None
    if center:
        xy_prime = R @ (xy - centroid)
    else:
        xy_prime = R @ (xy - centroid) + centroid

    return xy_prime

def lascsv_transform(incsvfname, outcsvfname, shpfname,
    groupby, xcolumn, ycolumn,
    theta, center = False, delimiter = ',', verbose = False):
    # load files
    in_df = pandas.read_csv(incsvfname, sep = delimiter)
    shpFile = shapefile.Reader(shpfname)

    # drop na from in_df
    in_df = in_df.dropna()

    # iterate through shape entries
    for shp in shpFile:
        # if the shape is not a polygon (code == 5), skip it.
        if shp.shape.shapeType != 5:
            print(
                "skipping %s '%s' (%s != 5)" %
                (shp.shape.shapeTypeName, shp.record.id, shp.shape.shapeType)
            )
            continue

        # get the shape ID (plot ID)
        shpID = shp.record.id

        # construct a polygon for the shapefile entry
        polygon = Polygon(shp.shape.points)

        # extract centroid for the polygon
        centroid = numpy.array(polygon.centroid.coords).T # (2,1)

        # get mask of where shpID == groupby column value
        mask = in_df.loc[:,groupby] == shpID

        # get x vector
        x_vec = in_df.loc[mask,xcolumn].values

        # get y vector
        y_vec = in_df.loc[mask,ycolumn].values

        # construct xy matrix
        xy_mat = numpy.stack([x_vec, y_vec])

        # calculate rotated matrix
        xy_prime_mat = rotate_mat(
            xy = xy_mat,
            centroid = centroid,
            theta = theta,
            center = center
        )

        # update x position coordinates
        in_df.loc[mask,xcolumn] = xy_prime_mat[0,:]

        # update y position coordinates
        in_df.loc[mask,ycolumn] = xy_prime_mat[1,:]

        if verbose:
            print("Processed %s" % shpID)

    # save to file
    in_df.to_csv(outcsvfname, sep = delimiter, index = None)



if __name__ == '__main__':
    # build parser object
    parser = optparse.OptionParser()

    # add options to parser
    parser.add_option(
        "-i", "--incsv",
        dest = "incsvfname",
        help = "Mandatory input CSV file name.",
        metavar = "FILE"
    )
    parser.add_option(
        "-o", "--outcsv",
        dest = "outcsvfname",
        help = "Mandatory output CSV file name.",
        metavar = "FILE"
    )
    parser.add_option(
        "-s", "--shapefile",
        dest = "shpfname",
        help = "Mandatory input SHP file name.",
        metavar = "FILE"
    )
    parser.add_option(
        "-x", "--xcolumn",
        dest = "xcolumn",
        help = "Column name of x values.",
        metavar = "STR"
    )
    parser.add_option(
        "-y", "--ycolumn",
        dest = "ycolumn",
        help = "Column name of y values.",
        metavar = "STR"
    )
    parser.add_option(
        "-g", "--groupby",
        dest = "groupby",
        help = "Column name to group by.",
        metavar = "STR"
    )
    parser.add_option(
        "-t", "--theta",
        dest = "theta",
        help = "Angle to rotate. This angle must be in radians!",
        metavar = "FLOAT"
    )
    parser.add_option(
        "-d", "--delimiter",
        dest = "delimiter",
        help = "Specify a custom delimiter for the CSV.",
        metavar = "STR"
    )
    parser.add_option(
        "-c", "--center",
        dest = "center",
        help = "Center polygon centroid coordinates to the origin.",
        action = "store_true"
    )
    parser.add_option(
        "-v", "--verbose",
        dest = "verbose",
        help = "Process data verbosely.",
        action = "store_true"
    )

    # parse the arguments
    options, args = parser.parse_args()

    if not options.incsvfname:
        parser.error("Input CSV file name not given.")
    if not os.path.isfile(options.incsvfname):
        parser.error("Input CSV file does not exist.")
    if not options.outcsvfname:
        parser.error("Output CSV file name not given.")
    if not options.shpfname:
        parser.error("Input SHP file name not given.")
    if not os.path.isfile(options.shpfname):
        parser.error("Input SHP file does not exist.")
    if not options.groupby:
        parser.error("No name grouping column given.")
    if not options.xcolumn:
        parser.error("No x column name given.")
    if not options.ycolumn:
        parser.error("No y column name given.")
    if not options.theta:
        parser.error("No rotation angle provided.")
    if not options.delimiter:
        options.delimiter = ','

    # transform the data
    lascsv_transform(
        incsvfname = options.incsvfname,
        outcsvfname = options.outcsvfname,
        shpfname = options.shpfname,
        groupby = options.groupby,
        xcolumn = options.xcolumn,
        ycolumn = options.ycolumn,
        theta = float(options.theta),
        center = options.center,
        delimiter = options.delimiter,
        verbose = options.verbose
    )
