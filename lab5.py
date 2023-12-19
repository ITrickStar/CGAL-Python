import os
from CGAL.CGAL_Polyhedron_3 import Polyhedron_3
from PyQt5.QtCore import *
from PyQt5.QtWidgets import *
from PyQGLViewer import *
from OpenGL.GL import *
import argparse

from CGAL.CGAL_Polyhedron_3 import Polyhedron_3
from CGAL.CGAL_Mesh_3 import Polyhedral_mesh_domain_3
from CGAL.CGAL_Mesh_3 import Mesh_3_parameters
from CGAL.CGAL_Mesh_3 import Default_mesh_criteria
from CGAL import CGAL_Mesh_3


class Viewer(QGLViewer):
    def __init__(self, path, facet_size, parent=None):
        QGLViewer.__init__(self, parent)
        self.vertices = []  # : [Point_3()]
        self.polyhedron = Polyhedron_3()
        self.path = path
        self.facet_size = facet_size
        self.c3t3 = None

    def draw(self):
        # glPointSize(2.5)

        # glBegin(GL_POINTS)
        # glColor3f(0.0, 1.0 , 0.0)
        # for vertex in self.vertices:
        #     glVertex3f( vertex.x(),
        #                 vertex.y(),
        #                 vertex.z())
        # glEnd()
        self.finite_element_mesh_3d()

        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE)
        for j, facet in enumerate(self.polyhedron.facets()):
            glBegin(GL_POLYGON)
            glColor3f(1, 0.5, 0.5)
            circulator = facet.facet_begin()
            for i in range(facet.facet_degree()):
                v = circulator.next()
                glVertex3f(v.vertex().point().x(),
                           v.vertex().point().y(),
                           v.vertex().point().z())
            glEnd()

        if self.c3t3 is not None:
            glPointSize(2.5)

            for j, cell in enumerate(self.c3t3.cells()):
                glBegin(GL_POLYGON)
                if j == 0:
                    glColor3f(0.0, 1.0, 0.0)
                else:
                    glColor3f(0.0, 0.0, 1.0)
                for i in range(4):
                    for k in range(3):
                        glVertex3f(cell.vertex(i).point().x(),
                                   cell.vertex(i).point().y(),
                                   cell.vertex(i).point().z())
                        glVertex3f(cell.vertex(k).point().x(),
                                   cell.vertex(k).point().y(),
                                   cell.vertex(k).point().z())
                glEnd()
        else:
            print("c3t3 is None")

    def keyPressEvent(self, e):
        print("keyPressEvent")
        if (e.nativeVirtualKey() == Qt.Key_1):
            self.finite_element_mesh_3d()
        elif (e.nativeVirtualKey() == Qt.Key_C):
            self.clear()
        self.updateGL()

    def finite_element_mesh_3d(self):
        self.polyhedron = Polyhedron_3(self.path)
        polyhedrons = []
        polyhedrons.append(self.polyhedron)

        domain = Polyhedral_mesh_domain_3(*polyhedrons)
        params = Mesh_3_parameters()
        criteria = Default_mesh_criteria()
        criteria.facet_size(self.facet_size)

        self.c3t3 = CGAL_Mesh_3.make_mesh_3(domain, criteria, params)
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE)

    def clear(self):
        print("clear")
        self.c3t3 = None
        self.updateGL()


def main():
    parser = argparse.ArgumentParser(
        description='Draw a 3D mesh.')
    parser.add_argument('input', type=str, help='Input mesh file')
    parser.add_argument('-fs', '--facet_size', type=float, default=0.15,  # 0.15
                        help='This parameter controls the size of surface facets. Each surface facet has a surface '
                             'Delaunay ball which is a ball circumscribing the surface facet and centered on the '
                             'surface patch. The parameter facet_size is either a constant or a spatially variable '
                             'scalar field, providing an upper bound for the radii of surface Delaunay balls. Set to -1 to disable.')
    args = parser.parse_args()

    obj_to_draw = args.input
    facet_size = args.facet_size
    if not os.path.exists(obj_to_draw):
        raise Exception('Input file \''+input+'\' does not exist.')

    qapp = QApplication([])
    viewer = Viewer(path=obj_to_draw, facet_size=facet_size)
    viewer.show()
    qapp.exec_()


if __name__ == '__main__':
    main()
