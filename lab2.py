# convex_hull_2d.py
from __future__ import print_function
from CGAL.CGAL_Kernel import Point_2, Point_3
from CGAL import CGAL_Convex_hull_2, CGAL_Convex_hull_3
import random
import numpy as np
from PyQt5.QtCore import *
from PyQt5.QtWidgets import *
from PyQGLViewer import *
from OpenGL.GL import *
from OpenGL.GLUT import *
from CGAL.CGAL_Polyhedron_3 import Polyhedron_3


def print_2d_points(points):
    for p in points:
        print(f"({p})", end=" ")
        print()


class MyObject:
    def __init__(self, vertices, color=(1, 0, 1)):
        self.vertices = vertices
        self.color = color


class Viewer(QGLViewer):
    figures = []

    def __init__(self, *args, parent=None):
        QGLViewer.__init__(self, parent)
        self.figures.extend(*args)

    # def initializeGL(self):
    #     glEnable(GL_DEPTH_TEST)

    def draw_2D(self, obj):
        glBegin(GL_POINTS)

        result = []
        CGAL_Convex_hull_2.convex_hull_2(obj.vertices, result)
        inversed = tuple(1 - x for x in obj.color)
        glColor3f(*inversed)
        for vertex in result:
            x, y = vertex.x(), vertex.y()
            glVertex2f(x, y)

        glColor3f(*obj.color)
        for vertex in obj.vertices:
            x, y = vertex.x(), vertex.y()
            glVertex2f(x, y)
        glEnd()

        glBegin(GL_LINE_LOOP)
        glColor3f(*inversed)
        for vertex in result:
            x, y = vertex.x(), vertex.y()
            glVertex2f(x, y)
        glEnd()

    def draw_3D(self, obj):
        glBegin(GL_POINTS)

        result = Polyhedron_3()
        CGAL_Convex_hull_3.convex_hull_3(obj.vertices, result)
        inversed = tuple(1 - x for x in obj.color)
        glColor3f(*inversed)
        for vertex in result.facets():
            facet = vertex.facet_begin()
            for v in range(vertex.facet_degree()):
                v = facet.next()
                x, y, z = (
                    v.vertex().point().x(),
                    v.vertex().point().y(),
                    v.vertex().point().z(),
                )
                glVertex3f(x, y, z)

        glColor3f(*obj.color)
        for vertex in obj.vertices:
            x, y, z = vertex.x(), vertex.y(), vertex.z()
            glVertex3f(x, y, z)
        glEnd()

        glColor3f(*inversed)
        for vertex in result.facets():
            glBegin(GL_LINE_LOOP)
            facet = vertex.facet_begin()
            for v in range(vertex.facet_degree()):
                v = facet.next()
                x, y, z = (
                    v.vertex().point().x(),
                    v.vertex().point().y(),
                    v.vertex().point().z(),
                )
                glVertex3f(x, y, z)
            glEnd()

    def draw(self):
        glPointSize(5.0)
        for obj in self.figures:
            if hasattr(obj.vertices[0], "z"):
                self.draw_3D(obj)
            else:
                self.draw_2D(obj)

    def clear_vertices(self):
        self.figures.clear()
        self.updateGL()

    def keyPressEvent(self, e):
        modifiers = e.modifiers()
        if e.nativeVirtualKey() == Qt.Key_W:
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE)
        elif e.nativeVirtualKey() == Qt.Key_F:
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL)
        elif e.nativeVirtualKey() == Qt.Key_C:
            self.clear_vertices()

        self.updateGL()


def main():
    L2 = []
    points = 250
    for i in range(points):
        L2.append(Point_2(random.random(), random.random()))

    L3 = []
    for i in range(points):
        L3.append(Point_3(random.random(), random.random(), random.random()))

    figs = []
    obj1 = MyObject(vertices=L2, color=(1, 0, 1))
    figs.append(obj1)
    obj2 = MyObject(vertices=L3, color=(1, 0, 1))
    figs.append(obj2)

    qapp = QApplication([])
    viewer = Viewer(figs)
    viewer.show()
    qapp.exec_()


if __name__ == "__main__":
    main()
