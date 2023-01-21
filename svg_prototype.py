#!/usr/bin/env python3

"""

1. Input an SVG
2. Convert the SVG to a polygon
3. Find a bounding box for the polygon
4. For the studs, iterate over where the center of each stud would be and calculate the shortest distance
   to the polygon. If it's less than the stud radius, then we know a stud can't fit
   If it's less than the stud radius plus the wall width plus our tolerances, it should fit if we stack the piece
   directly on top with 
5. If the "notches" are requested on underside, we need to do the same as #4 but cut out 
6. We will need to generate negative space using some method of generating an inset
7. For each post, iterate over where teh center of each post would be and calculate the shortest distance to the
    polygon. If it's less than the post radius, we know a post can't fit
"""

from pathlib import Path

# Read SVG into a list of path objects and list of dictionaries of attributes
# from svgpathtools import svg2paths, wsvg
import click
from pydantic import BaseModel

from svgpath2mpl import parse_path
import numpy as np

from shapely import Polygon, bounds, Point
from shapely import affinity

import math
import solid


class BlockDimensions(BaseModel):
    # Attributes could be Decimal
    post_wall_thickness: float
    wall_thickness: float
    stud_diameter: float
    hollow_stud_inner_diameter: float
    stud_height: float
    stud_spacing: float

    unit_height: float
    baseplate_height: float

    pin_diameter: float
    post_diameter: float
    cylinder_precision: float
    reinforcing_width: float

    spline_length: float
    spline_thickness: float

    horizontal_hole_diameter: float
    horizontal_hole_z_offset: float
    horizontal_hole_bevel_diameter: float
    horizontal_hole_bevel_depth: float

    axle_spline_width: float
    axle_diameter: float

    wall_play: float = 0.1
    horizontal_hole_wall_thickness = 1.0


lego_block_dimensions = BlockDimensions(
    post_wall_thickness=0.85,
    wall_thickness=1.45,
    stud_diameter=4.85,  # Radius is 2.425
    hollow_stud_inner_diameter=3.1,
    stud_height=1.8,
    stud_spacing=8.0,
    unit_height=9.6,
    baseplate_height=1.3,
    pin_diameter=3.0,
    post_diameter=6.5,
    cylinder_precision=0.1,
    reinforcing_width=0.7,
    spline_length=0.25,
    spline_thickness=0.7,
    horizontal_hole_diameter=4.8,
    horizontal_hole_z_offset=5.8,
    horizontal_hole_bevel_diameter=6.2,
    horizontal_hole_bevel_depth=0.9,
    axle_spline_width=2.0,
    axle_diameter=5,
)

# How close can a stud be to an edge
# Not sure this is right
STUD_EDGE_CLEARANCE = 0.075


def generate_regular_polygon(n_sides):
    theta = np.linspace(0, 2 * np.pi, n_sides, endpoint=False)

    # This rotates the polygon so that the top and bottom are parallel if there are an even number of sides
    theta += np.pi / n_sides

    return Polygon(np.stack([np.sin(theta), np.cos(theta)], axis=1))


def _render_polygon(polygon: Polygon):
    rings = [polygon.exterior, *polygon.interiors]

    points = []
    paths = []

    for ring in rings:
        path = []

        for x, y in ring.coords:
            path.append(len(points))
            points.append([x, y])

        paths.append(path)

    return solid.polygon(points, paths=paths)


def generate_swept_back_wing(*, block_dimensions: BlockDimensions):
    polygon = Polygon([[0, 0], [1, 0], [11, 1], [12, 1], [12, 7], [10, 7], [0, 2]])

    return affinity.scale(
        polygon, block_dimensions.stud_spacing, block_dimensions.stud_spacing
    )


def generate_rectangle(width=2, length=4, *, block_dimensions):
    polygon_coordinates = [
        (0, 0),
        (0, width),
        (length, width),
        (length, 0),
    ]

    return solid.polygon(polygon_coordinates)


def generate_block(
    polygon: Polygon,
    height: float = 1,
    stud_type: str = "solid",
    block_bottom_type: str = "open",
    include_wall_splines: bool = True,
    horizontal_holes: bool = False,
    vertical_axle_holes: bool = False,
    reinforcement: bool = False,  # TODO Make 'None' with default
    stud_notches: bool = False,
    round_radius: float = 0,
    stud_rescale: float = 1,
    stud_top_roundness: float = 0,
    dual_sided: bool = False,
    dual_bottom: bool = False,
    *,
    block_dimensions
):
    pass


ACID_PATH = """
    m 19.054693,25.397583 v 10.683594 h 6.328125 q 3.183594,0 4.707031,-1.308594 1.542969,-1.328125 1.542969,-4.042969 0,-2.734375 -1.542969,-4.023437 -1.523437,-1.308594 -4.707031,-1.308594 z m 0,-11.992187 v 8.789062 h 5.839844 q 2.890625,0 4.296875,-1.074219 1.425781,-1.09375 1.425781,-3.320312 0,-2.207031 -1.425781,-3.300781 -1.40625,-1.09375 -4.296875,-1.09375 z m -3.945312,-3.242188 h 10.078125 q 4.511718,0 6.953125,1.875 2.441406,1.875 2.441406,5.332031 0,2.675782 -1.25,4.257813 -1.25,1.582031 -3.671875,1.972656 2.910156,0.625 4.511719,2.617188 1.621093,1.972656 1.621093,4.941406 0,3.90625 -2.65625,6.035156 -2.65625,2.128906 -7.558593,2.128906 h -10.46875 z
"""

SWEPT_WING_PATH = """
    M 74.4693 159.512 L 74.4693 148.215 L 74.9426 146.771 L 75.697 145.728 L 76.6209 144.993 L 77.6026 144.469 L 159.436 107.482 L 166.51 107.482 L 166.51 152.057 L 159.436 152.057 L 133.272 154.681 L 89.6092 158.882 L 81.9868 159.507 L 74.4693 159.512
"""


def _svg_path_to_polygon(path: str) -> Polygon:

    parsed = parse_path(path)

    # Is it guaranteed that holes and outline are in this order
    *holes, outline = parsed.to_polygons()

    polygon = Polygon(outline, holes=holes)

    # I think SVG's coordinate system and openscad consider X coordinates backwards from one another
    polygon = affinity.scale(polygon, -1, 1)

    return polygon


# Scales polygon so that its dimensions are ready to turn into a block
def _normalize_polygon(
    polygon: Polygon, *, input_dimensions: str = "mm", block_dimensions: BlockDimensions
) -> Polygon:
    assert input_dimensions in {"mm", "studs"}

    scale = 1 if input_dimensions == "studs" else block_dimensions.stud_spacing

    assert input_dimensions == "mm"

    min_x, min_y, max_x, max_y = polygon.bounds

    width = max_x - min_x
    desired_width = scale * round(width / scale)
    length = max_y - min_y
    desired_length = scale * round(length / scale)

    polygon = affinity.translate(polygon, -min_x, -min_y)
    polygon = affinity.scale(
        polygon, desired_width / width, desired_length / length, origin=(0, 0)
    )

    return polygon


def _get_inner_points(polygon, *, block_dimensions: BlockDimensions):
    min_x, min_y, max_x, max_y = polygon.bounds

    rings = [polygon.exterior, *polygon.interiors]

    def edge_distance(point):
        return min(ring.distance(point) for ring in rings)

    stud_radius = block_dimensions.stud_diameter / 2
    required_clearance = stud_radius + STUD_EDGE_CLEARANCE

    for x in np.arange(
        block_dimensions.stud_spacing / 2,
        max_x,
        block_dimensions.stud_spacing,
    ):
        for y in np.arange(
            block_dimensions.stud_spacing / 2,
            max_y,
            block_dimensions.stud_spacing,
        ):
            point = Point(x, y)

            if polygon.contains(point) and edge_distance(point) >= required_clearance:
                yield x, y


def _get_overlapping_grid_points(polygon, *, block_dimensions: BlockDimensions):
    min_x, min_y, max_x, max_y = polygon.bounds

    rings = [polygon.exterior, *polygon.interiors]

    def edge_distance(point):
        return min(ring.distance(point) for ring in rings)

    stud_radius = block_dimensions.stud_diameter / 2
    required_clearance = (
        stud_radius + STUD_EDGE_CLEARANCE + block_dimensions.wall_thickness
    )

    for x in np.arange(
        -block_dimensions.stud_spacing / 2,
        max_x + block_dimensions.stud_spacing,
        block_dimensions.stud_spacing,
    ):
        for y in np.arange(
            -block_dimensions.stud_spacing / 2,
            max_y + block_dimensions.stud_spacing,
            block_dimensions.stud_spacing,
        ):
            point = Point(x, y)

            if edge_distance(point) < required_clearance:
                yield x, y


def _get_post_locations(polygon, *, block_dimensions: BlockDimensions):
    min_x, min_y, max_x, max_y = polygon.bounds

    rings = [polygon.exterior, *polygon.interiors]

    def edge_distance(point):
        return min(ring.distance(point) for ring in rings)

    post_radius = block_dimensions.post_diameter / 2
    required_clearance = (
        post_radius + STUD_EDGE_CLEARANCE + block_dimensions.post_wall_thickness
    )

    for x in np.arange(
        0,
        max_x + block_dimensions.stud_spacing / 2,  # Can be any epsilon
        block_dimensions.stud_spacing,
    ):
        for y in np.arange(
            0,
            max_y + block_dimensions.stud_spacing / 2,  # Can be any epsilon
            block_dimensions.stud_spacing,
        ):
            point = Point(x, y)

            if polygon.contains(point) and edge_distance(point) >= required_clearance:
                yield x, y


def inset(polygon, amount):
    buffer = 0.1

    min_x, min_y, max_x, max_y = polygon.bounds

    obj = _render_polygon(polygon)

    mold = (
        solid.translate([min_x - buffer, min_y - buffer])(
            solid.square(
                [2 * buffer + max_x - min_x, 2 * buffer + max_y - min_y], center=False
            )
        )
        - obj
    )

    return solid.intersection()([solid.minkowski()(mold, solid.circle(r=amount)), obj])


def _render_brick(polygon, *, block_dimensions: BlockDimensions):

    height_in_block_units = 1 / 3  # Make a param
    height = height_in_block_units * block_dimensions.unit_height

    stud = solid.cylinder(
        r=block_dimensions.stud_diameter / 2,
        h=block_dimensions.stud_height,
        segments=60,  # TODO Optimize this
    )

    polygon = _normalize_polygon(polygon, block_dimensions=block_dimensions)

    NOTCHES = True

    notches = []

    if NOTCHES:

        notch = solid.circle(
            r=block_dimensions.stud_diameter / 2,
            segments=60,  # TODO Optimize this
        )

        # We could do a 'crosses boundary' check here if we wanted to make our generate scad more terse
        for x, y in _get_overlapping_grid_points(
            polygon, block_dimensions=block_dimensions
        ):
            notches.append(solid.translate([x, y])(notch))

    walls = solid.linear_extrude(height=height)(
        inset(polygon, block_dimensions.wall_thickness)
    )

    # Strange place for the constant
    roof_thickness = 1

    roof = solid.translate([0, 0, height - roof_thickness])(
        solid.linear_extrude(height=roof_thickness)(_render_polygon(polygon))
    )

    shell = roof + walls

    studs = []

    for x, y in _get_inner_points(polygon, block_dimensions=block_dimensions):
        studs.append(solid.translate([x, y])(stud))

    if studs:
        shell += solid.translate([0, 0, height])(sum(studs))

    # Notches, make optional
    if notches:
        shell -= sum(notches)

    # TODO Add precision to cylinders
    # TODO Add axle
    post = solid.linear_extrude(height=height)(
        solid.circle(
            # No $fs support?
            d=block_dimensions.post_diameter,
            segments=60,
        )
        - solid.circle(
            d=block_dimensions.post_diameter - 2 * block_dimensions.post_wall_thickness,
            segments=60,
        )
    )

    for x, y in _get_post_locations(polygon, block_dimensions=block_dimensions):
        shell += solid.translate([x, y])(post)

    return shell


def _swept_wing_from_svg():
    polygon = _svg_path_to_polygon(SWEPT_WING_PATH)

    return polygon


@click.command()
@click.argument("output_file", type=Path)
def main(output_file):
    # Approximate circle of radius 2 block units
    # polygon = generate_regular_polygon(8)
    # polygon = affinity.scale(
    #     polygon,
    #     2 * lego_block_dimensions.stud_spacing,
    #     2 * lego_block_dimensions.stud_spacing,
    # )

    # polygon = generate_swept_back_wing(block_dimensions=lego_block_dimensions)

    polygon = _swept_wing_from_svg()

    with output_file.open("wt") as f:
        f.write(
            solid.scad_render(
                _render_brick(polygon, block_dimensions=lego_block_dimensions)
            ),
        )

    return

    assert False, polygon

    with open("rectangle_lego.scad", "wt") as f:
        root = generate_rectangle(block_dimensions=lego_block_dimensions)

        f.write(solid.scad_render(root))

    return

    studs = [
        solid.translate([x, y])(solid.cylinder(r=1, h=1))
        for x, y in _iterate_points(polygon)
    ]

    assert studs

    with open("tmp.scad", "wt") as f:
        f.write(solid.scad_render(SVGImport("tmp.svg", center=True) + sum(studs)))


if __name__ == "__main__":
    main()
