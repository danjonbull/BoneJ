package org.doube.bonej.geomorph;

import java.awt.Color;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.util.ArrayList;

import javax.media.j3d.Transform3D;
import javax.media.j3d.TransformGroup;
import javax.media.j3d.View;
import javax.vecmath.Color3f;
import javax.vecmath.Point3d;
import javax.vecmath.Point3f;
import javax.vecmath.Vector3d;

import customnode.CustomLineMesh;

import ij.IJ;
import ij.gui.Roi;
import ij3d.Content;
import ij3d.Image3DUniverse;
import ij3d.UniverseListener;

/**
 * 
 * 
 * @author Michael Doube
 */
public class Crosshairs implements UniverseListener, KeyListener {
	private static final Transform3D emptyTransform = new Transform3D();
	private Image3DUniverse univ;
	private double x, y, z;
	private CustomLineMesh clmX, clmY, clmZ;
	private Content cX, cY, cZ;
	private TransformGroup tg = new TransformGroup();
	private Transform3D t1 = new Transform3D();
	private Vector3d vector = new Vector3d();
	private Point3d globalMax, globalMin;

	public Crosshairs(double x, double y, double z, Image3DUniverse univ) {
		this.x = x;
		this.y = y;
		this.z = z;
		this.univ = univ;
		create();
	}

	private void create() {
		globalMax = new Point3d();
		univ.getGlobalMaxPoint(globalMax);

		globalMin = new Point3d();
		univ.getGlobalMinPoint(globalMin);

		ArrayList<Point3f> meshZ = new ArrayList<Point3f>();
		Point3f start1 = new Point3f();
		start1.x = (float) x;
		start1.y = (float) y;
		start1.z = (float) globalMin.z;
		meshZ.add(start1);

		Point3f end1 = new Point3f();
		end1.x = (float) x;
		end1.y = (float) y;
		end1.z = (float) globalMax.z;
		meshZ.add(end1);

		ArrayList<Point3f> meshY = new ArrayList<Point3f>();
		Point3f start2 = new Point3f();
		start2.x = (float) x;
		start2.y = (float) globalMin.y;
		start2.z = (float) z;
		meshY.add(start2);

		Point3f end2 = new Point3f();
		end2.x = (float) x;
		end2.y = (float) globalMax.y;
		end2.z = (float) z;
		meshY.add(end2);

		ArrayList<Point3f> meshX = new ArrayList<Point3f>();
		Point3f start3 = new Point3f();
		start3.x = (float) globalMin.x;
		start3.y = (float) y;
		start3.z = (float) z;
		meshX.add(start3);

		Point3f end3 = new Point3f();
		end3.x = (float) globalMax.x;
		end3.y = (float) y;
		end3.z = (float) z;
		meshX.add(end3);

		clmX = new customnode.CustomLineMesh(meshX);
		clmY = new customnode.CustomLineMesh(meshY);
		clmZ = new customnode.CustomLineMesh(meshZ);
		try {
			cX = univ.addCustomMesh(clmX, "Crosshairs X");
			cY = univ.addCustomMesh(clmY, "Crosshairs Y");
			cZ = univ.addCustomMesh(clmZ, "Crosshairs Z");
			resetColor();
			cX.setLocked(true);
			cY.setLocked(true);
			cZ.setLocked(true);
		} catch (NullPointerException npe) {
			IJ.log("3D Viewer was closed before rendering completed.");
			return;
		} catch (Exception e) {
			IJ.log(e.getMessage());
		}
		univ.addUniverseListener(this);
		univ.getCanvas().addKeyListener(this);
	}

	public void set(double x, double y, double z) {
		this.x = x;
		this.y = y;
		this.z = z;
	}

	public Point3d get() {
		return new Point3d(x, y, z);
	}

	public void update() {

		if (x < globalMin.x)
			x = globalMin.x;
		if (x > globalMax.x)
			x = globalMax.x;
		if (y < globalMin.y)
			y = globalMin.y;
		if (y > globalMax.y)
			y = globalMax.y;
		if (z < globalMin.z)
			z = globalMin.z;
		if (z > globalMax.z)
			z = globalMax.z;

		Point3f start1 = new Point3f();
		start1.x = (float) x;
		start1.y = (float) y;
		start1.z = (float) globalMin.z;

		Point3f end1 = new Point3f();
		end1.x = (float) x;
		end1.y = (float) y;
		end1.z = (float) globalMax.z;

		Point3f start2 = new Point3f();
		start2.x = (float) x;
		start2.y = (float) globalMin.y;
		start2.z = (float) z;

		Point3f end2 = new Point3f();
		end2.x = (float) x;
		end2.y = (float) globalMax.y;
		end2.z = (float) z;

		Point3f start3 = new Point3f();
		start3.x = (float) globalMin.x;
		start3.y = (float) y;
		start3.z = (float) z;

		Point3f end3 = new Point3f();
		end3.x = (float) globalMax.x;
		end3.y = (float) y;
		end3.z = (float) z;

		clmZ.setCoordinate(0, start1);
		clmZ.setCoordinate(1, end1);
		clmY.setCoordinate(0, start2);
		clmY.setCoordinate(1, end2);
		clmX.setCoordinate(0, start3);
		clmX.setCoordinate(1, end3);

		cX.setTransform(emptyTransform);
		cY.setTransform(emptyTransform);
		cZ.setTransform(emptyTransform);

	}

	public void show() {
		cX.setVisible(true);
		cY.setVisible(true);
		cZ.setVisible(true);
	}

	public void hide() {
		cX.setVisible(false);
		cY.setVisible(false);
		cZ.setVisible(false);
	}

	private void resetColor() {
		Color c = Roi.getColor();
		Color3f color = new Color3f((float)(c.getRed()/255),
				(float)(c.getGreen()/255), (float)(c.getBlue()/255));
		cX.setColor(color);
		cY.setColor(color);
		cZ.setColor(color);
	}

	
	public void canvasResized() {
		// TODO Auto-generated method stub

	}

	
	public void contentAdded(Content c) {
		// TODO Auto-generated method stub

	}

	
	public void contentChanged(Content c) {
	}

	
	public void contentRemoved(Content c) {
		// TODO Auto-generated method stub

	}

	
	public void contentSelected(Content c) {
		resetColor();
		Color3f red = new Color3f(1.0f, 0.0f, 0.0f);
		try {
			if (c.equals(cX)) {
				cX.setColor(red);
				univ.setShowBoundingBoxUponSelection(false);
			} else if (c.equals(cY)) {
				cY.setColor(red);
				univ.setShowBoundingBoxUponSelection(false);
			} else if (c.equals(cZ)) {
				cZ.setColor(red);
				univ.setShowBoundingBoxUponSelection(false);
			} else {
				univ.setShowBoundingBoxUponSelection(true);
				return;
			}
		} catch (Exception e) {
			IJ.log("Please update your 3D Viewer to" +
					" at least 22.09.2010 version");
			IJ.log(e.getMessage());
		}
	}

	
	public void transformationFinished(View view) {
		update();
	}

	
	public void transformationStarted(View view) {
		univ.getGlobalMaxPoint(globalMax);
		univ.getGlobalMinPoint(globalMin);
	}

	
	public void transformationUpdated(View view) {
		try {
			if (univ.getSelected().equals(cX)) {
				tg = cX.getLocalTranslate();
				tg.getTransform(t1);
				t1.get(vector);
				y += vector.y;
				z += vector.z;
			} else if (univ.getSelected().equals(cY)) {
				tg = cY.getLocalTranslate();
				tg.getTransform(t1);
				t1.get(vector);
				x += vector.x;
				z += vector.z;
			} else if (univ.getSelected().equals(cZ)) {
				tg = cZ.getLocalTranslate();
				tg.getTransform(t1);
				t1.get(vector);
				x += vector.x;
				y += vector.y;
			}
			update();
		} catch (Exception e) {
			IJ.log(e.getMessage());
		}

	}

	
	public void universeClosed() {
		// TODO Auto-generated method stub

	}

	
	public void keyPressed(KeyEvent e) {
		// allow translation when shift is pressed
		if (e.getKeyCode() == 16) {
			cX.setLocked(false);
			cY.setLocked(false);
			cZ.setLocked(false);
		}
	}

	
	public void keyReleased(KeyEvent e) {
		// prevent rotation
		if (e.getKeyCode() == 16) {
			resetColor();
			cX.setSelected(false);
			cY.setSelected(false);
			cZ.setSelected(false);
			cX.setLocked(true);
			cY.setLocked(true);
			cZ.setLocked(true);
		}
	}

	
	public void keyTyped(KeyEvent e) {
		// TODO Auto-generated method stub

	}
}
