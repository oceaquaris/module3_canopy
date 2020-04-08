#!/usr/bin/env python3
import optparse
import numpy
import pandas
import shapefile
import os
from shapely.geometry.polygon import Polygon
from scipy.interpolate import bisplrep, bisplev

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

def lascsv_transform(in_df, shpFile, groupby, xcolumn, ycolumn, theta,
    center = False, verbose = False):
    # drop na from in_df
    in_df = in_df.dropna()

    # make a dictionary of centroid
    centroid_dict = {
        None : numpy.array([[0.0],[0.0]])
    }

    # iterate through shape entries, getting entries and their centroids
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

        # add key to centroid_dict
        centroid_dict[shpID] = centroid

        if verbose:
            print("Calculated centroid %s" % shpID)

    # construct centroid_list
    centroid_list = []
    for group in in_df.loc[:,groupby].values:
        centroid_list.append(centroid_dict.get(group, centroid_dict[None]))

    # get (2,p) matrix of centroid
    centroid_mat = numpy.concatenate(centroid_list, axis = 1)

    # get x vector
    x_vec = in_df.loc[:,xcolumn].values

    # get y vector
    y_vec = in_df.loc[:,ycolumn].values

    # construct xy matrix
    xy_mat = numpy.stack([x_vec, y_vec])

    # calculate rotated matrix
    xy_prime_mat = rotate_mat(
        xy = xy_mat,
        centroid = centroid_mat,
        theta = theta,
        center = center
    )

    # update x position coordinates
    # raises warning and suggests me to use what I am currently doing. wtf
    in_df.loc[:,xcolumn] = xy_prime_mat[0,:]

    # update y position coordinates
    in_df.loc[:,ycolumn] = xy_prime_mat[1,:]

    if verbose:
        print("Transformed positions")

    print(in_df)

    # return the transformed file
    return in_df

def lascsv_interpolate(in_df, shpFile, groupby, xcolumn, ycolumn,
    xstart, xstop, xnum, ystart, ystop, ynum, kind = "linear", verbose = False):
    # build components
    xcomp = numpy.linspace(xstart, xstop, xnum)
    ycomp = numpy.linspace(ystart, ystop, ynum)

    # build meshgrid
    mesh = numpy.meshgrid(xcomp, ycomp, indexing = 'xy')

    # extract x,y mesh
    xmesh = mesh[0].flatten()
    ymesh = mesh[1].flatten()

    df_list = []

    # shape index
    shpIX = 0

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

        # get mask of where shpID == groupby column value
        mask = in_df.loc[:,groupby] == shpID

        # get x vector
        x_vec = in_df.loc[mask,xcolumn].values

        # get y vector
        y_vec = in_df.loc[mask,ycolumn].values

        # make dictionary
        df_dict = {
            "x_position": xmesh,
            "y_position": ymesh
        }

        # calculate z values separately
        tck = bisplrep(x_vec, y_vec, in_df.loc[mask,"z_position"], kx=1, ky=1, s=40)
        df_dict["z_position"] = bisplev(xcomp, ycomp, tck).flatten()

        # for values to interpolate:
        for col in ["r_record", "g_record", "b_record"]:
            # build model
            tck = bisplrep(x_vec, y_vec, in_df.loc[mask,col], kx=1, ky=1, s=40)

            # interpolate and add to dictionary
            df_dict[col] = bisplev(xcomp, ycomp, tck).flatten()

        # add shape ID and IX to dictionary
        df_dict["shpID"] = shpID
        df_dict["shpIX"] = shpIX

        # make pandas.DataFrame and add to list
        df_list.append(pandas.DataFrame(df_dict))

        # increment shape index
        shpIX += 1

        if verbose:
            print("Interpolated %s" % shpID)

    # concatenate all DataFrame together
    out_df = pandas.concat(df_list)

    return out_df

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
        "--xstart",
        dest = "xstart",
        help = "Starting x location to interpolate.",
        metavar = "FLOAT"
    )
    parser.add_option(
        "--xstop",
        dest = "xstop",
        help = "Stopping x location to interpolate.",
        metavar = "FLOAT"
    )
    parser.add_option(
        "--xnum",
        dest = "xnum",
        help = "Number of x locations to interpolate.",
        metavar = "FLOAT"
    )
    parser.add_option(
        "-y", "--ycolumn",
        dest = "ycolumn",
        help = "Column name of y values.",
        metavar = "STR"
    )
    parser.add_option(
        "--ystart",
        dest = "ystart",
        help = "Starting y location to interpolate.",
        metavar = "FLOAT"
    )
    parser.add_option(
        "--ystop",
        dest = "ystop",
        help = "Stopping y location to interpolate.",
        metavar = "FLOAT"
    )
    parser.add_option(
        "--ynum",
        dest = "ynum",
        help = "Number of y locations to interpolate.",
        metavar = "FLOAT"
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
        "-k", "--kind",
        dest = "kind",
        help = "Interpolation kind",
        metavar = "STR"
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
    if not options.xstart:
        parser.error("No x start position given.")
    if not options.xstop:
        parser.error("No x stop position given.")
    if not options.xnum:
        parser.error("No x number of locations given.")
    if not options.ycolumn:
        parser.error("No y column name given.")
    if not options.ystart:
        parser.error("No y start position given.")
    if not options.ystop:
        parser.error("No y stop position given.")
    if not options.ynum:
        parser.error("No y number of locations given.")
    if not options.theta:
        parser.error("No rotation angle provided.")
    if not options.kind:
        options.kind = "linear"
    if not options.delimiter:
        options.delimiter = ','

    # load files
    df = pandas.read_csv(options.incsvfname, sep = options.delimiter)
    shpFile = shapefile.Reader(options.shpfname)

    # transform the data
    df = lascsv_transform(
        in_df = df,
        shpFile = shpFile,
        groupby = options.groupby,
        xcolumn = options.xcolumn,
        ycolumn = options.ycolumn,
        theta = float(options.theta),
        center = options.center,
        verbose = options.verbose
    )

    # interpolate the data
    df = lascsv_interpolate(
        in_df = df,
        shpFile = shpFile,
        groupby = options.groupby,
        xcolumn = options.xcolumn,
        ycolumn = options.ycolumn,
        xstart = float(options.xstart),
        xstop = float(options.xstop),
        xnum = int(options.xnum),
        ystart = float(options.ystart),
        ystop = float(options.ystop),
        ynum = int(options.ynum),
        kind = options.kind,
        verbose = options.verbose
    )

    # write to file
    df.to_csv(outcsvfname, sep = delimiter, index = None)
