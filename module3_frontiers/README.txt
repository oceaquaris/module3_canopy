These files are CSV files *with a header*.

The header is mostly self-explanatory as to what is in the column.

shpIX : an integer index of the polygon grouping.
shpID : an associated string with that polygon (i.e. plot name)

A couple of notes:
    Some of the points have a shpIX and shpID of 4294967295 and NA, respectively.
        This indicates that the point is not assigned to a specific plot polygon.
        Remark: This may be a result of an approximation ray-tracing algorithm in the prefiltration step.
    RGB values are not in the range [0,255] because they are multiplied by 256.
        This is for compatibility with cameras that can record more than 256 RGB values.
        To get RGB, divide by 256
    XYZ values are in ft. For Z this is feet above sea level.
