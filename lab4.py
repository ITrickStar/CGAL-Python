#!/usr/bin/python3

import os
import sys
import numpy as np
from tempfile import NamedTemporaryFile
import meshio
from CGAL.CGAL_Polyhedron_3 import Polyhedron_3
from CGAL.CGAL_Mesh_3 import Mesh_3_Complex_3_in_triangulation_3
from CGAL.CGAL_Mesh_3 import Polyhedral_mesh_domain_3
from CGAL.CGAL_Mesh_3 import Mesh_3_parameters
from CGAL.CGAL_Mesh_3 import Default_mesh_criteria
from CGAL import CGAL_Mesh_3


def print_stats(c3t3, m):
    triangulation = c3t3.triangulation()
    print(F"  Number of points     : {triangulation.number_of_vertices():,}".replace(
        ',', ' '))
    print(F"  Number of facets     : {triangulation.number_of_facets():,}".replace(
        ',', ' '))
    print(F"  Number of cells      : {triangulation.number_of_cells():,}".replace(
        ',', ' '))

    v = 0.
    n_degenerate = 0
    for c in c3t3.cells():
        t = triangulation.tetrahedron(c)
        v += t.volume()
        if t.is_degenerate():
            n_degenerate += 1
    print(
        F"  Total volume         : {v / 1000:.2f} ml (approx. {v * 1.07 / 1000 / 1000:.2f} kg)")
    print(F"  Number of degenerate : {n_degenerate}")

    # stats
    points = m.points
    tetras = None
    if isinstance(m.cells, list):
        for cell in m.cells:
            if cell.type == 'tetra':
                tetras = cell.data
    else:
        if 'tetra' in m.cells:
            tetras = m.cells['tetra']

    # https://www2.eecs.berkeley.edu/Pubs/TechRpts/2007/EECS-2007-104.pdf
    angles = np.zeros(len(tetras))
    for i in range(len(tetras)):
        t = tetras[i]
        nodes = points[t]
        u = nodes[1] - nodes[0]
        v = nodes[2] - nodes[0]
        w = nodes[3] - nodes[0]
        cos_theta = (np.dot(np.cross(u, v), np.cross(u, w))) / \
            (np.linalg.norm(np.cross(u, v)) * np.linalg.norm(np.cross(u, w)))
        angles[i] = cos_theta
    angles = np.rad2deg(np.arccos(angles))
    print(F"  Mean dihedral angle  : {np.mean(angles):.2f}")
    print(F"  Median dihedral angle: {np.median(angles):.2f}")
    print(F"  Min dihedral angle   : {np.min(angles):.2f}")
    print(F"  Max dihedral angle   : {np.max(angles):.2f}")
    print(F"  Std dihedral angle   : {np.std(angles):.2f}")


def triangulate(
        input, input_meshes, verbose,
        facet_angle, facet_size, facet_distance,
        cell_radius_edge_ratio, cell_size,
        lloyd_iterations, lloyd_time, odt_iterations, odt_time,
        perturber_time, perturber_sliver_bound, exudation_time, exudation_sliver_bound):

    # Load input meshes
    output = "D:\\Works\\CGAL\\Python_labs\\log.off"
    polyhedrons = []
    # for input_mesh in input_meshes:
    # with NamedTemporaryFile() as tempfile:
    # meshio.write(tempfile, input_mesh, file_format="off")

    # Create input polyhedron
    # polyhedron = Polyhedron_3(tempfile.name)
    meshi o.write(output, *input_meshes, file_format="off")
    polyhedron = Polyhedron_3(input)
    polyhedrons.append(polyhedron)

    if verbose:
        for i, polyhedron in zip(range(len(polyhedrons)), polyhedrons):
            print(F"Input mesh {i+1} has:")
            print(
                F"  Number of points : {polyhedron.size_of_vertices():,}".replace(',', ' '))
            print(
                F"  Number of facets : {polyhedron.size_of_facets():,}".replace(',', ' '))

    # Create domain
    if verbose:
        print("Creating domain...", end='', flush=True)
    domain = Polyhedral_mesh_domain_3(*polyhedrons)
    if verbose:
        print(" Done", flush=True)

    params = Mesh_3_parameters()
    if odt_iterations > -1:
        params.set_odt(odt_time, odt_iterations, 0.02, 0.01)
    if lloyd_iterations > -1:
        params.set_lloyd(lloyd_time, lloyd_iterations, 0.02, 0.01)
    if perturber_time > -1:
        params.set_perturb(perturber_time, perturber_sliver_bound)
    if exudation_time > -1:
        params.set_exude(exudation_time, exudation_sliver_bound)

    # Mesh criteria (no cell_size set)
    criteria = Default_mesh_criteria()
    if facet_angle >= 0:
        criteria.facet_angle(facet_angle)
    if facet_distance >= 0:
        criteria.facet_distance(facet_distance)
    if facet_size >= 0:
        criteria.facet_size(facet_size)
    if cell_radius_edge_ratio >= 0:
        criteria.cell_radius_edge_ratio(cell_radius_edge_ratio)
    if cell_size >= 0:
        criteria.cell_size(cell_size)

    # Mesh generation
    verbose and print("Generating the mesh...", end='', flush=True)
    c3t3 = CGAL_Mesh_3.make_mesh_3(domain, criteria, params)
    verbose and print(" Done", flush=True)

    output_mesh = None
    # with NamedTemporaryFile() as tempfile:
    
    c3t3.output_to_medit(output)
    output_mesh = meshio.read(output, file_format="medit")
    # output_mesh = meshio.read(output, file_format="off")

    if verbose:
        print(F"Statistics on output mesh:")
        print_stats(c3t3, output_mesh)

    return output_mesh

    # meshio.write(output_filename, m)


def __main():
    import argparse

    parser = argparse.ArgumentParser(
        description='Triangulate a 3D surface mesh.')
    parser.add_argument('input', type=str, nargs='+', help='Input mesh file')
    parser.add_argument('output', type=str, help='Output mesh file')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Output process informations')
    parser.add_argument('-fa', '--facet_angle', type=float, default=-1,  # 25
                        help='This parameter controls the shape of surface facets. Specifically, it is a lower bound '
                             'for the angle (in degrees) of surface facets. When boundary surfaces are smooth, the '
                             'termination of the meshing process is guaranteed if this angular bound is at most 30 '
                             'degrees. Set to -1 to disable.')
    parser.add_argument('-fs', '--facet_size', type=float, default=-1,  # 0.15
                        help='This parameter controls the size of surface facets. Each surface facet has a surface '
                             'Delaunay ball which is a ball circumscribing the surface facet and centered on the '
                             'surface patch. The parameter facet_size is either a constant or a spatially variable '
                             'scalar field, providing an upper bound for the radii of surface Delaunay balls. Set to -1 to disable.')
    parser.add_argument('-fd', '--facet_distance', type=float, default=-1,  # 0.008
                        help='This parameter controls the approximation error of boundary and subdivision surfaces. '
                             'Specifically, it is either a constant or a spatially variable scalar field. It provides '
                             'an upper bound for the distance between the circumcenter of a surface facet and the '
                             'center of a surface Delaunay ball of this facet. Set to -1 to disable.')
    parser.add_argument('-cr', '--cell_radius_edge_ratio', type=float, default=-1,  # 3
                        help='This parameter controls the shape of mesh cells (but can\'t filter slivers). It is an '
                             'upper bound for the ratio between the circumradius of a mesh tetrahedron and its '
                             'shortest edge. There is a theoretical bound for this parameter: the Delaunay refinement p'
                             'rocess is guaranteed to terminate for values of cell_radius_edge_ratio bigger than 2. Set to -1 to disable.')
    parser.add_argument('-cs', '--cell_size', type=float, default=-1,  # 1
                        help='This parameter controls the size of mesh tetrahedra. It is either a scalar or a spatially '
                             'variable scalar field. It provides an upper bound on the circumradii of the mesh '
                             'tetrahedra. Set to -1 to disable.')
    parser.add_argument('-li', '--lloyd_iterations', type=int, default=0,
                        help="Maximum number of iteration for the global Lloyd optimization pass. Set to -1 to disable. Set to 0 for no limit.")
    parser.add_argument('-lt', '--lloyd_time', type=int, default=0,
                        help="Maximum number of time in seconds for the global Lloyd optimization pass. Set to 0 for no limit.")
    parser.add_argument('-oi', '--odt_iterations', type=int, default=0,
                        help="Maximum number of iteration for the global ODT-smoother optimization pass. Set to -1 to disable. Set to 0 for no limit.")
    parser.add_argument('-ot', '--odt_time', type=int, default=0,
                        help="Maximum number of time in seconds for the global ODT-smoother optimization pass. Set to 0 for no limit.")
    parser.add_argument('-pt', '--perturber_time', type=int, default=0,
                        help="Maximum number of time in seconds for the local perturber optimization pass. Set to -1 to disable. Set to 0 for no limit.")
    parser.add_argument('-ps', '--perturber_sliver_bound', type=float, default=0,
                        help="Targeted lower bound in degrees on dihedral angles of mesh cells. The optimization will "
                             "runs as long as steps are successful and step number sliver_bound (after which the worst "
                             "tetrahedron in the mesh has a smallest angle larger than sliver_bound degrees) has not "
                             "been reached. The default value is 0 and means that there is no targeted bound: the "
                             "perturber then runs as long as steps are successful.")
    parser.add_argument('-et', '--exudation_time', type=int, default=0,
                        help="Maximum number of time in seconds for the local exudation process. Set to -1 to disable. Set to 0 for no limit.")
    parser.add_argument('-es', '--exudation_sliver_bound', type=float, default=0,
                        help="Targeted lower bound in degrees on dihedral angles of mesh cells. The exudation process "
                             "considers in turn all the mesh cells that have a smallest dihedral angle less than "
                             "sliver_bound and tries to make them disappear by weighting their vertices. "
                             "The optimization process stops when every cell in the mesh achieves this quality. "
                             "The default value is 0 and means that there is no targeted bound: the exuder then runs as "
                             "long as it can improve the smallest dihedral angles of the set of cells incident to some vertices.")
    args = parser.parse_args()

    input_meshes = []
    for input in args.input:
        if not os.path.exists(input):
            raise Exception('Input file \''+input+'\' does not exist.')
        input_meshes.append(meshio.read(input))

    output_mesh = triangulate(input, input_meshes, args.verbose, args.facet_angle, args.facet_size, args.facet_distance,
                              args.cell_radius_edge_ratio, args.cell_size,
                              args.lloyd_iterations, args.lloyd_time, args.odt_iterations, args.odt_time,
                              args.perturber_time, args.perturber_sliver_bound, args.exudation_time, args.exudation_sliver_bound)

    # meshio.write(args.output, output_mesh)
    meshio.write(args.output, output_mesh, file_format="medit")

if __name__ == "__main__":
    __main()
