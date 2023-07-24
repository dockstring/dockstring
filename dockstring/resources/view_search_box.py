# Inspired by drawBoundingBox.py by JasonVertrees
# PyMOL script for visualization of bounding box

from pymol import cmd
from pymol.cgo import BEGIN, COLOR, END, LINES, LINEWIDTH, VERTEX


def view_search_box(center_x, center_y, center_z, size_x, size_y, size_z):
    line_width = 2.0
    rgb = (1.0, 1.0, 1.0)

    # Convert center_{x,y,z} and size_{x,y,z} to floats
    center_x = float(center_x)
    center_y = float(center_y)
    center_z = float(center_z)
    size_x = float(size_x)
    size_y = float(size_y)
    size_z = float(size_z)

    # Compute min and max coordinates of bound box
    min_x = center_x - size_x / 2
    min_y = center_y - size_y / 2
    min_z = center_z - size_z / 2
    max_x = center_x + size_x / 2
    max_y = center_y + size_y / 2
    max_z = center_z + size_z / 2

    # Define the search box
    # yapf: disable
    search_box = [
        LINEWIDTH, line_width,

        BEGIN, LINES,
        COLOR, *rgb,

        VERTEX, min_x, min_y, min_z,  # 1
        VERTEX, min_x, min_y, max_z,  # 2

        VERTEX, min_x, max_y, min_z,  # 3
        VERTEX, min_x, max_y, max_z,  # 4

        VERTEX, max_x, min_y, min_z,  # 5
        VERTEX, max_x, min_y, max_z,  # 6

        VERTEX, max_x, max_y, min_z,  # 7
        VERTEX, max_x, max_y, max_z,  # 8

        VERTEX, min_x, min_y, min_z,  # 1
        VERTEX, max_x, min_y, min_z,  # 5

        VERTEX, min_x, max_y, min_z,  # 3
        VERTEX, max_x, max_y, min_z,  # 7

        VERTEX, min_x, max_y, max_z,  # 4
        VERTEX, max_x, max_y, max_z,  # 8

        VERTEX, min_x, min_y, max_z,  # 2
        VERTEX, max_x, min_y, max_z,  # 6

        VERTEX, min_x, min_y, min_z,  # 1
        VERTEX, min_x, max_y, min_z,  # 3

        VERTEX, max_x, min_y, min_z,  # 5
        VERTEX, max_x, max_y, min_z,  # 7

        VERTEX, min_x, min_y, max_z,  # 2
        VERTEX, min_x, max_y, max_z,  # 4

        VERTEX, max_x, min_y, max_z,  # 6
        VERTEX, max_x, max_y, max_z,  # 8

        END
    ]
    # yapf: enable

    cmd.load_cgo(search_box, 'search_box')
    return 'search_box'


cmd.extend('view_search_box', view_search_box)
