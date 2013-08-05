package org.bonej;

/**
 * EllipsoidFactor plugin for ImageJ
 * Copyright 2013 Michael Doube
 * 
 *This program is free software: you can redistribute it and/or modify
 *it under the terms of the GNU General Public License as published by
 *the Free Software Foundation, either version 3 of the License, or
 *(at your option) any later version.
 *
 *This program is distributed in the hope that it will be useful,
 *but WITHOUT ANY WARRANTY; without even the implied warranty of
 *MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *GNU General Public License for more details.
 *
 *You should have received a copy of the GNU General Public License
 *along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

import javax.vecmath.Color3f;
import javax.vecmath.Point3f;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.plugin.PlugIn;
import ij.process.ByteProcessor;
import ij.gui.GenericDialog;
import ij.macro.Interpreter;
import ij.measure.Calibration;
import ij3d.Image3DUniverse;

import org.doube.geometry.Vectors;
import org.doube.geometry.Ellipsoid;
import org.doube.jama.EigenvalueDecomposition;
import org.doube.jama.Matrix;
import org.doube.skeleton.Skeletonize3D;
import org.doube.util.ArrayHelper;
import org.doube.util.ImageCheck;
import org.doube.util.ResultInserter;
import org.doube.util.UsageReporter;

import customnode.CustomPointMesh;

/**
 * <p>
 * <b>Plate_Rod</b>
 * </p>
 * <p>
 * ImageJ plugin to describe the local geometry of a binary image in an
 * oblate/prolate spheroid space. Uses Skeletonize3D to generate a 3D skeleton,
 * the points of which are used as centres for star volumes. Local geometry is
 * determined by the ratio between the first and second eigenvalues and first
 * and third eigenvalues of each star volume.
 * </p>
 * 
 * @author Michael Doube
 * 
 */
public class EllipsoidFactor implements PlugIn, Comparator<Ellipsoid> {
	private final byte foreground = (byte) 255;
	private int nVectors = 100;
	private Image3DUniverse universe = new Image3DUniverse();

	/**
	 * increment for vector searching in real units. Defaults to ~Nyquist
	 * sampling of a unit pixel
	 */
	private double vectorIncrement = 1 / 2.3;

	/**
	 * Number of skeleton points per ellipsoid. Sets the granularity of the
	 * ellipsoid fields.
	 */
	private int skipRatio = 50;

	public void run(String arg) {
		if (!ImageCheck.checkEnvironment())
			return;
		ImagePlus imp = IJ.getImage();
		if (imp == null) {
			IJ.noImage();
			return;
		}
		ImageCheck ic = new ImageCheck();
		if (!ic.isBinary(imp) || !ic.isMultiSlice(imp)
				|| !ic.isVoxelIsotropic(imp, 0.001)) {
			IJ.error("8-bit binary stack with isotropic pixel spacing required.");
			return;
		}
		Calibration cal = imp.getCalibration();
		String units = cal.getUnits();
		vectorIncrement *= Math.min(cal.pixelDepth,
				Math.min(cal.pixelHeight, cal.pixelWidth));
		GenericDialog gd = new GenericDialog("Setup");
		gd.addNumericField("Sampling increment", vectorIncrement, 3, 8, units);
		gd.addNumericField("Vectors", nVectors, 0, 8, "");
		gd.addNumericField("Skeleton points per ellipsoid", skipRatio, 0);
		gd.addHelp("http://bonej.org/ef");
		gd.showDialog();
		if (!Interpreter.isBatchMode()) {
			vectorIncrement = gd.getNextNumber();
			nVectors = (int) Math.round(gd.getNextNumber());
			skipRatio = (int) Math.round(gd.getNextNumber());
		}
		if (gd.wasCanceled())
			return;

		final double[][] unitVectors = Vectors.regularVectors(nVectors);
		int[][] skeletonPoints = skeletonPoints(imp);

		IJ.log("Found " + skeletonPoints.length + " skeleton points");

		universe.show();

		Ellipsoid[] ellipsoids = findEllipsoids(imp, skeletonPoints,
				unitVectors);

		IJ.log("Found " + ellipsoids.length + " ellipsoids");

		float[][] biggestEllipsoid = findBiggestEllipsoid(imp, ellipsoids);

		ImageStack bigStack = new ImageStack(imp.getWidth(), imp.getHeight());
		for (int i = 1; i < biggestEllipsoid.length; i++)
			bigStack.addSlice("" + i, biggestEllipsoid[i]);

		ImagePlus bigImp = new ImagePlus("", bigStack);
		bigImp.setDisplayRange(-ellipsoids.length / 2, ellipsoids.length);
		bigImp.show();

		ResultInserter ri = ResultInserter.getInstance();
		ri.updateTable();
		UsageReporter.reportEvent(this).send();
	}

	/**
	 * For each foreground pixel of the input image, find the ellipsoid of
	 * greatest volume
	 * 
	 * @param imp
	 * @param ellipsoids
	 * @return array containing the indexes of the biggest ellipsoids which
	 *         contain each point
	 */
	private float[][] findBiggestEllipsoid(ImagePlus imp, Ellipsoid[] ellipsoids) {

		ImageStack stack = imp.getImageStack();
		final int w = stack.getWidth();
		final int h = stack.getHeight();
		final int d = stack.getSize();

		Calibration cal = imp.getCalibration();
		final double vW = cal.pixelWidth;
		final double vH = cal.pixelHeight;
		final double vD = cal.pixelDepth;

		float[][] biggest = new float[d + 1][w * h];

		for (int z = 1; z <= d; z++) {
			byte[] slicePixels = (byte[]) stack.getPixels(z);
			float[] bigSlice = biggest[z];
			// -1 means background, 0 will be the biggest ellipsoid
			Arrays.fill(bigSlice, -ellipsoids.length);
			for (int y = 0; y < h; y++) {
				int offset = y * w;
				for (int x = 0; x < w; x++) {
					if (slicePixels[offset + x] == foreground) {
						bigSlice[offset + x] = biggestEllipsoid(ellipsoids, x
								* vW, y * vH, z * vD);
					}
				}
			}
		}

		return biggest;
	}

	/**
	 * Search the list of ellipsoids and return the index of the largest
	 * ellipsoid which contains the point x, y, z
	 * 
	 * @param ellipsoids
	 * @param x
	 * @param y
	 * @param z
	 * @return the index of the largest ellipsoid which contains this point
	 */
	private int biggestEllipsoid(Ellipsoid[] ellipsoids, double x, double y,
			double z) {
		final int l = ellipsoids.length;
		for (int i = 0; i < l; i++) {
			if (ellipsoids[i].contains(x, y, z))
				return i;
		}
		return -1;
	}

	/**
	 * Using skeleton points as seeds, propagate along each vector until a
	 * boundary is hit. Use the resulting cloud of boundary points as input into
	 * an ellipsoid fit.
	 * 
	 * @param imp
	 * @param skeletonPoints
	 * @param unitVectors
	 * @return
	 */
	private Ellipsoid[] findEllipsoids(ImagePlus imp, int[][] skeletonPoints,
			double[][] unitVectors) {
		final int nPoints = skeletonPoints.length;
		Ellipsoid[] ellipsoids = new Ellipsoid[nPoints];

		// make sure array contains null in the non-calculated elements
		Arrays.fill(ellipsoids, null);

		for (int i = 0; i < nPoints; i += skipRatio) {
			IJ.showStatus("Optimising ellipsoid " + (i + 1) + "/" + nPoints);
			ellipsoids[i] = optimiseEllipsoid(imp, skeletonPoints[i],
					unitVectors);
		}

		ellipsoids = ArrayHelper.removeNulls(ellipsoids);

		// Sort using this class' compare method
		Arrays.sort(ellipsoids, this);
		return ellipsoids;
	}

	/**
	 * given a seed point, find the ellipsoid which best fits the binarised
	 * structure
	 * 
	 * @param imp
	 * @param is
	 * @param unitVectors
	 * @return ellipsoid fitting the point cloud of boundaries lying at the end
	 *         of vectors surrounding the seed point. If ellipsoid fitting
	 *         fails, returns null
	 */
	private Ellipsoid optimiseEllipsoid(final ImagePlus imp,
			int[] skeletonPoint, double[][] unitVectors) {

		Calibration cal = imp.getCalibration();
		final double pW = cal.pixelWidth;
		final double pH = cal.pixelHeight;
		final double pD = cal.pixelDepth;

		ImageStack stack = imp.getImageStack();

		// cache slices into an array
		ByteProcessor[] ips = new ByteProcessor[stack.getSize() + 1];
		for (int i = 1; i <= stack.getSize(); i++) {
			ips[i] = (ByteProcessor) stack.getProcessor(i);
		}

		final int w = ips[1].getWidth();
		final int h = ips[1].getHeight();
		final int d = ips.length - 1;

		// centre point of vector field
		final double px = skeletonPoint[0] * pW;
		final double py = skeletonPoint[1] * pH;
		final double pz = skeletonPoint[2] * pD;

		// Instantiate a small spherical ellipsoid
		final double[][] orthogonalVectors = { { 1, 0, 0 }, { 0, 1, 0 },
				{ 0, 0, 1 } };

		Ellipsoid ellipsoid = new Ellipsoid(0, 0, 0, px, py, pz,
				orthogonalVectors);

		// dilate the sphere until it hits the background
		while (isContained(ellipsoid, ips, pW, pH, pD, w, h, d)) {
			ellipsoid.dilate(vectorIncrement);
			IJ.showStatus("Ellipsoid volume = " + ellipsoid.getVolume());
		}

		IJ.log("Sphere fit with radius " + ellipsoid.getMajorRadius());

		// get the points of contact
		List<double[]> contactPoints = findSphereContactPoints(ellipsoid, ips,
				pW, pH, pD, w, h, d);

		// calculate the new orientation matrix for the ellipsoid
		// short axis aligned to the contact points
		ellipsoid = reOrientEllipsoid(contactPoints, ellipsoid);

		// display axes in 3D viewer
		drawAxes(ellipsoid);

		// add points of contact to the 3D viewer
		List<Point3f> contactPointsf = new ArrayList<Point3f>(
				contactPoints.size());
		for (double[] p : contactPoints) {
			Point3f point = new Point3f((float) p[0], (float) p[1],
					(float) p[2]);
			contactPointsf.add(point);
		}

		// contract the ellipsoid by one increment so all points are
		// inside the foregrounds
		// ellipsoid.contract(vectorIncrement);

		double[][] pointCloud = ellipsoid.getSurfacePoints(100);

		// List<Point3f> pointList = new ArrayList<Point3f>();
		// for (int p = 0; p < pointCloud.length; p++) {
		// if (pointCloud[p] == null)
		// continue;
		// Point3f e = new Point3f();
		// e.x = (float) pointCloud[p][0];
		// e.y = (float) pointCloud[p][1];
		// e.z = (float) pointCloud[p][2];
		// pointList.add(e);
		// }
		// CustomPointMesh mesh = new CustomPointMesh(pointList);
		// mesh.setPointSize(2.0f);
		Color3f cColour = new Color3f((float) (px / pW) / w, (float) (py / pH)
				/ h, (float) (pz / pD) / d);
		// mesh.setColor(cColour);

		CustomPointMesh contactPointMesh = new CustomPointMesh(contactPointsf);
		contactPointMesh.setPointSize(2.5f);
		Color3f invColour = new Color3f(1 - cColour.x, 1 - cColour.y,
				1 - cColour.z);
		contactPointMesh.setColor(invColour);

		try {
			// universe.addCustomMesh(mesh,
			// "Point cloud " + px + " " + py + " " + pz).setLocked(true);
			universe.addCustomMesh(contactPointMesh,
					"Contact points of " + px + " " + py + " " + pz).setLocked(
					true);

		} catch (NullPointerException npe) {
			IJ.log("3D Viewer was closed before rendering completed.");
		}

		// now dilate long and middle axes until they hit the sides
		// new points will be added to contact point list

		return ellipsoid;
	}

	private void drawAxes(Ellipsoid ellipsoid) {
		final double[] centroid = ellipsoid.getCentre();
		final double cX = centroid[0];
		final double cY = centroid[1];
		final double cZ = centroid[2];

		final double[] radii = ellipsoid.getRadii();
		final double ra = radii[0];
		final double rb = radii[1];
		final double rc = radii[2];

		final double[][] eV = ellipsoid.getEigenVectors();

		List<Point3f> longAxis = new ArrayList<Point3f>();
		List<Point3f> middleAxis = new ArrayList<Point3f>();
		List<Point3f> shortAxis = new ArrayList<Point3f>();

		Point3f start1 = new Point3f();
		start1.x = (float) (cX - ra * eV[0][0]);
		start1.y = (float) (cY - ra * eV[1][0]);
		start1.z = (float) (cZ - ra * eV[2][0]);
		longAxis.add(start1);

		Point3f end1 = new Point3f();
		end1.x = (float) (cX + ra * eV[0][0]);
		end1.y = (float) (cY + ra * eV[1][0]);
		end1.z = (float) (cZ + ra * eV[2][0]);
		longAxis.add(end1);

		Point3f start2 = new Point3f();
		start2.x = (float) (cX - rb * eV[0][1]);
		start2.y = (float) (cY - rb * eV[1][1]);
		start2.z = (float) (cZ - rb * eV[2][1]);
		middleAxis.add(start2);

		Point3f end2 = new Point3f();
		end2.x = (float) (cX + rb * eV[0][1]);
		end2.y = (float) (cY + rb * eV[1][1]);
		end2.z = (float) (cZ + rb * eV[2][1]);
		middleAxis.add(end2);

		Point3f start3 = new Point3f();
		start3.x = (float) (cX - rc * eV[0][2]);
		start3.y = (float) (cY - rc * eV[1][2]);
		start3.z = (float) (cZ - rc * eV[2][2]);
		shortAxis.add(start3);

		Point3f end3 = new Point3f();
		end3.x = (float) (cX + rc * eV[0][2]);
		end3.y = (float) (cY + rc * eV[1][2]);
		end3.z = (float) (cZ + rc * eV[2][2]);
		shortAxis.add(end3);

		Color3f red = new Color3f(1.0f, 0.0f, 0.0f);
		Color3f green = new Color3f(0.0f, 1.0f, 0.0f);
		Color3f blue = new Color3f(0.0f, 0.0f, 1.0f);
		try {
			universe.addLineMesh(longAxis, red, "Long axis " + cX + cY + cZ,
					false).setLocked(true);
			universe.addLineMesh(middleAxis, green,
					"Middle axis " + cX + cY + cZ, false).setLocked(true);
			universe.addLineMesh(shortAxis, blue, "Short axis " + cX + cY + cZ,
					false).setLocked(true);
		} catch (NullPointerException npe) {
			IJ.log("3D Viewer was closed before rendering completed.");
			return;
		}

	}

	private Ellipsoid reOrientEllipsoid(List<double[]> contactPoints,
			Ellipsoid ellipsoid) {
		// generate a covariance matrix from the contact points
		final double[] centroid = ellipsoid.getCentre();
		final double cX = centroid[0];
		final double cY = centroid[1];
		final double cZ = centroid[2];

		double sDxDx = 0;
		double sDyDy = 0;
		double sDzDz = 0;
		double sDxDy = 0;
		double sDxDz = 0;
		double sDyDz = 0;
		for (double[] p : contactPoints) {
			final double dx = p[0] - cX;
			final double dy = p[1] - cY;
			final double dz = p[2] - cZ;
			sDxDx += dx * dx;
			sDyDy += dy * dy;
			sDzDz += dz * dz;
			sDxDy += dx * dy;
			sDxDz += dx * dz;
			sDyDz += dy * dz;
		}
		double[][] c = new double[3][3];
		c[0][0] = sDxDx;
		c[1][1] = sDyDy;
		c[2][2] = sDzDz;
		c[0][1] = sDxDy;
		c[0][2] = sDxDz;
		c[1][0] = sDxDy;
		c[1][2] = sDyDz;
		c[2][0] = sDxDz;
		c[2][1] = sDyDz;
		// do the eigendecomposition, maybe SVD?

		Matrix C = new Matrix(c);
		System.out.println("Covariance matrix");
		C.print(10, 8);

		EigenvalueDecomposition E = C.eig();

		System.out.println("Eigenvectors");
		E.getV().print(10, 8);

		System.out.println("Eigenvalues");
		E.getD().print(10, 8);

		// assign the largest eigenvalue's vector to be the short axis of the
		// ellipsoid

		ellipsoid.setEigenVectors(E.getV().getArrayCopy());

		return ellipsoid;
	}

	private List<double[]> findSphereContactPoints(Ellipsoid ellipsoid,
			ByteProcessor[] ips, final double pW, final double pH,
			final double pD, final int w, final int h, final int d) {
		double[][] points = ellipsoid.getSurfacePoints(nVectors);

		List<double[]> contactPoints = new ArrayList<double[]>();

		for (double[] p : points) {
			final int x = (int) Math.floor(p[0] / pW);
			final int y = (int) Math.floor(p[1] / pH);
			final int z = (int) Math.floor(p[2] / pD);
			if (isOutOfBounds(x, y, z, w, h, d))
				continue;
			if ((byte) ips[z].get(x, y) != foreground)
				contactPoints.add(p);
		}
		return contactPoints;
	}

	private boolean isContained(Ellipsoid ellipsoid, ByteProcessor[] ips,
			final double pW, final double pH, final double pD, final int w,
			final int h, final int d) {
		double[][] points = ellipsoid.getSurfacePoints(nVectors);
		for (double[] p : points) {
			final int x = (int) Math.floor(p[0] / pW);
			final int y = (int) Math.floor(p[1] / pH);
			final int z = (int) Math.floor(p[2] / pD);
			if (isOutOfBounds(x, y, z, w, h, d))
				continue;
			if ((byte) ips[z].get(x, y) != foreground)
				return false;
		}
		return true;
	}

	/**
	 * return true if pixel coordinate is out of image bounds
	 * 
	 * @param x
	 * @param y
	 * @param z
	 * @param w
	 * @param h
	 * @param d
	 * @return
	 */
	private boolean isOutOfBounds(int x, int y, int z, int w, int h, int d) {
		if (x < 0 || x >= w || y < 0 || y >= h || z < 1 || z > d)
			return true;
		else
			return false;
	}

	private int[][] skeletonPoints(ImagePlus imp) {
		Skeletonize3D sk = new Skeletonize3D();
		ImagePlus skeleton = sk.getSkeleton(imp);
		ImageStack skeletonStack = skeleton.getStack();

		skeleton.show();

		final int d = imp.getStackSize();
		final int h = imp.getHeight();
		final int w = imp.getWidth();

		IJ.log("Skeleton image is " + w + " x " + h + " x " + d);

		ArrayList<int[]> list = new ArrayList<int[]>();

		for (int z = 1; z <= d; z++) {
			byte[] slicePixels = (byte[]) skeletonStack.getPixels(z);
			for (int y = 0; y < h; y++) {
				int offset = y * w;
				for (int x = 0; x < w; x++) {
					if (slicePixels[offset + x] == foreground) {
						int[] array = { x, y, z };
						list.add(array);
					}
				}
			}
		}

		IJ.log("Skeleton point ArrayList contains " + list.size() + " points");

		int[][] skeletonPoints = list.toArray(new int[list.size()][]);

		return skeletonPoints;
	}

	/**
	 * Compare Ellipsoids by volume.
	 * 
	 * Sorting based on this method will result in Ellipsoids sorted in order of
	 * <b>descending</b> volume.
	 * 
	 */
	public int compare(Ellipsoid o1, Ellipsoid o2) {
		return Double.compare(o2.getVolume(), o1.getVolume());
	}

}
