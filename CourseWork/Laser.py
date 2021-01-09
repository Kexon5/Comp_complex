import vector3d.vector as vec
import math
import numpy as np


def createPlane(line1, line2, point1, isLine=True):
    if isLine:
        point2 = point1 + line1
        point3 = point1 + line2
    else:
        point2 = line1
        point3 = line2

    det_x = (point2.y - point1.y) * (point3.z - point1.z) - (point3.y - point1.y) * (point2.z - point1.z)
    det_y = (point2.x - point1.x) * (point3.z - point1.z) - (point3.x - point1.x) * (point2.z - point1.z)
    det_z = (point2.x - point1.x) * (point3.y - point1.y) - (point3.x - point1.x) * (point2.y - point1.y)

    return [det_x, -det_y, det_z, -det_x * point1.x + det_y * point1.y - det_z * point1.z]


def checkInside(pointCheck, point1, point2, pointStart, isDetector=False):
    vecCheck = vec.Vector(pointCheck[0], pointCheck[1], pointCheck[2])
    vecCheck = (vecCheck - pointStart).normalize()
    if isDetector:
        vec1 = point1
        vec2 = point2
    else:
        vec1 = (point1 - pointStart).normalize()
        vec2 = (point2 - pointStart).normalize()

    x_solve = None
    flag = False
    b = np.array([vecCheck.x, vecCheck.y])
    A = np.array([[vec1.x, vec2.x], [vec1.y, vec2.y]])
    try:
        x_solve = np.linalg.solve(A, b)
    except np.linalg.LinAlgError:
        flag = True

    if flag:
        flag = False
        b = np.array([vecCheck.y, vecCheck.z])
        A = np.array([[vec1.y, vec2.y], [vec1.z, vec2.z]])
        try:
            x_solve = np.linalg.solve(A, b)
        except np.linalg.LinAlgError:
            flag = True

    if x_solve[0] * x_solve[1] >= 0:
        return True
    else:
        return False


def crossLinePlane(line, plane, lines_length=1e+2, isDetector=False):
    if isDetector:
        v = (line[1] - line[0]).normalize()
        numerator = plane[0] * line[1].x + plane[1] * line[1].y + plane[2] * line[1].z + plane[3]
        denominator = plane[0] * v.x + plane[1] * v.y + plane[2] * v.z
    else:
        numerator = plane[0] * line[1].x + plane[1] * line[1].y + plane[2] * line[1].z + plane[3]
        denominator = plane[0] * line[0].x + plane[1] * line[0].y + plane[2] * line[0].z
    if abs(denominator) <= 1e-5:
        return []
    else:
        t = -numerator / denominator
        if isDetector:
            vector = vec.Vector(line[1].x + v.x * t, line[1].y + v.y * t, line[1].z + v.z * t)
        else:
            vector = vec.Vector(line[1].x + line[0].x * t, line[1].y + line[0].y * t, line[1].z + line[0].z * t)
        if t < 0 or (line[1] - vector).length() > lines_length:
            return []
        else:
            return [vector.x, vector.y, vector.z]


class Laser:

    def __init__(self, center, normal, divergence_angle=0, lines_length=0):
        self.center = vec.Vector(center[0], center[1], center[2])
        self.normal = vec.Vector(normal[0], normal[1], normal[2]).normalize()

        v = vec.Vector(1.0, 0.0, 0.0)
        self.lines = [
            (v + vec.Vector(0.0, math.tan(divergence_angle), math.tan(divergence_angle))).normalize(),
            (v + vec.Vector(0.0, -math.tan(divergence_angle),
                            math.tan(divergence_angle))).normalize(),
            (v + vec.Vector(0.0, -math.tan(divergence_angle),
                            -math.tan(divergence_angle))).normalize(),
            (v + vec.Vector(0.0, math.tan(divergence_angle),
                            -math.tan(divergence_angle))).normalize()]

        def unit_vector(vector):
            return vector / np.linalg.norm(vector)

        def angle_between(v1, v2):
            v1_u = unit_vector(v1)
            v2_u = unit_vector(v2)
            return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

        rotVector = createPlane(v, self.normal, self.center)

        v_normal = np.array([self.normal.x, self.normal.y, self.normal.z])
        v_1 = np.array([v.x, v.y, v.z])

        angle = angle_between(v_1, v_normal)
        rotate_matrix = np.array([[math.cos(angle) + (1 - math.cos(angle)) * rotVector[0] ** 2,
                                   (1 - math.cos(angle)) * rotVector[0] * rotVector[1] - math.sin(angle) * rotVector[2],
                                   (1 - math.cos(angle)) * rotVector[0] * rotVector[2] + math.sin(angle) * rotVector[1]],
                                  [(1 - math.cos(angle)) * rotVector[0] * rotVector[1] + math.sin(angle) * rotVector[2],
                                   math.cos(angle) + (1 - math.cos(angle)) * rotVector[1] ** 2,
                                   (1 - math.cos(angle)) * rotVector[1] * rotVector[2] - math.sin(angle) * rotVector[0]],
                                  [(1 - math.cos(angle)) * rotVector[0] * rotVector[2] - math.sin(angle) * rotVector[1],
                                   (1 - math.cos(angle)) * rotVector[1] * rotVector[2] + math.sin(angle) * rotVector[0],
                                   math.cos(angle) + (1 - math.cos(angle)) * rotVector[2] ** 2]])
        for i in range(len(self.lines)):
            vector = np.array([self.lines[i].x, self.lines[i].y, self.lines[i].z])
            new_vec = np.dot(rotate_matrix, vector)
            self.lines[i] = vec.Vector(new_vec[0], new_vec[1], new_vec[2])

        self.lines_length = lines_length

    def plot(self, ax, **kwargs):
        for line in self.lines:
            p2 = self.center + line * self.lines_length
            ax.plot([self.center.x, p2.x], [self.center.y, p2.y], [self.center.z, p2.z], **kwargs)

    def intersection(self, detector, rows, cols, aperture):
        indexes = [0, rows - 1, rows * cols - 1, rows * (cols - 1)]
        ap = vec.Vector(aperture[0], aperture[1], aperture[2])
        point_list = [vec.Vector(detector.pixels[indexes[0]].pos[0], detector.pixels[indexes[0]].pos[1],
                                 detector.pixels[indexes[0]].pos[2]),
                      vec.Vector(detector.pixels[indexes[1]].pos[0], detector.pixels[indexes[1]].pos[1],
                                 detector.pixels[indexes[1]].pos[2]),
                      vec.Vector(detector.pixels[indexes[2]].pos[0], detector.pixels[indexes[2]].pos[1],
                                 detector.pixels[indexes[2]].pos[2]),
                      vec.Vector(detector.pixels[indexes[3]].pos[0], detector.pixels[indexes[3]].pos[1],
                                 detector.pixels[indexes[3]].pos[2])]

        detector_plane_list = []
        cross_list1 = []
        indexes = []
        for i in range(len(point_list)):
            detector_plane_list.append(
                createPlane(point_list[i], point_list[(i + 1) % len(point_list)], ap, isLine=False))
            indexes.clear()
            for j in range(len(self.lines)):

                elem = crossLinePlane([self.lines[j], self.center], detector_plane_list[i],
                                      lines_length=self.lines_length)

                if len(elem) == 3:
                    if checkInside(elem, point_list[i], point_list[(i + 1) % len(point_list)], ap):
                        cross_list1.append(elem)

        cross_list2 = []
        laser_plane = []
        for i in range(len(point_list)):
            laser_plane.append(createPlane(self.lines[i], self.lines[(i + 1) % len(point_list)], self.center))
            indexes.clear()
            for j in range(len(self.lines)):
                elem = crossLinePlane([point_list[j], ap], laser_plane[i],
                                      lines_length=self.lines_length, isDetector=True)

                if len(elem) == 3:
                    if checkInside(elem, self.lines[i], self.lines[(i + 1) % len(point_list)], self.center,
                                   isDetector=True):
                        cross_list2.append(elem)

        cross_list1 = cross_list1 + cross_list2
        cross_list1 = np.array(cross_list1)
        return cross_list1
