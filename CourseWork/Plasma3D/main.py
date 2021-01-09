import matplotlib.pyplot as plt
import numpy as np

from cylindrical_plasma import CylindricalPlasma, PlasmaDrawable
from detector import Detector, DetectorDrawable

from Computational_complexes.CourseWork.Laser import Laser

from scipy.spatial import ConvexHull, Delaunay


def example_plot():
    # creating plasma and 3 different detectors
    plasma = PlasmaDrawable(r_min=0.0, r_max=1.0, z_min=-1.0, z_max=1.0)
    plasma.build_segmentation(n_r=2, n_phi=4, n_z=2)
    d1 = DetectorDrawable(center=(-1.5, -1.5, 0.0), aperture=(-1.1, -1.0, 0.0), height=0.6, width=0.6)
    d1.set_pixels(rows=3, cols=3)
    l1 = Laser(center=(1.0, 0.0, 2.5), normal=(0.0, 0.0, -1.0), divergence_angle=np.radians(10), lines_length=1.1 * d1.lines_length_calculate(plasma))

    # setup plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.grid(False)
    ax.set_title('Example plot')

    points = l1.intersection(d1, d1.rows, d1.cols, d1.detector_metrics.aperture)
    d1.plot(ax, lines_length=d1.lines_length_calculate(plasma), color='green')
    l1.plot(ax, color='blue')

    hull_opt = ['Qr', 'Qs']

    for opt in hull_opt:
        hull = None
        try:
            hull = ConvexHull(points, qhull_options=opt)
        except:
            print()
        if hull != None:
            for simplex in hull.simplices:
                plt.plot(points[simplex, 0], points[simplex, 1], points[simplex, 2], 'r')
        else:
            for elem in points:
                ax.scatter(elem[0], elem[1], elem[2], color='red')

    plt.show()


example_plot()
