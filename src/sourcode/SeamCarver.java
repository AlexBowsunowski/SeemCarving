package sourcode;
import edu.princeton.cs.algs4.Picture;

import java.awt.Color;
import java.io.File;
import java.util.Arrays;

public class SeamCarver {

    private int[][] color;

    private double[][] energy;

    private double[][] distTo;
    private double distToSink;
    private int[][] edgeTo;
    private int edgeToSink;

    // size picture
    private int width;
    private int height;

    // create a seam carver object based on the given picture
    public SeamCarver(Picture picture) {
        if (picture == null) throw new java.lang.NullPointerException();

        width = picture.width();
        height = picture.height();

        color = new int[height][width];

        energy = new double[height][width];

        for (int i = 0; i < height; ++i) {
            for (int j = 0; j < width; ++j) {
                color[i][j] = picture.get(j, i).getRGB();
            }
        }

        for (int i = 0; i < height; ++i) {
            for (int j = 0; j < width; ++j) {
                energy[i][j] = calculateEnergy(j, i);
            }
        }
    }


    private double calculateEnergy(int x, int y) {
        if (x < 0 || y < 0 || x >= width() || y >= height())
            throw new IndexOutOfBoundsException();

        if (x == 0 || y == 0 || x == width() - 1 || y == height() - 1)
            return 1000;

        Color up = new Color(color[y - 1][x]);
        Color down = new Color(color[y + 1][x]);
        Color left = new Color(color[y][x - 1]);
        Color right = new Color(color[y][x + 1]);

        return Math.sqrt(gradient(up, down) + gradient(left, right));
    }

    private double gradient(Color a, Color b) {
        return Math.pow(a.getRed() - b.getRed(), 2) +
                Math.pow(a.getBlue() - b.getBlue(), 2) +
                Math.pow(a.getGreen() - b.getGreen(), 2);
    }

    // current picture
    public Picture picture() {
        Picture picture = new Picture(width(), height());
        for (int i = 0; i < height(); ++i)
            for (int j = 0; j < width(); ++j)
                picture.set(j, i, new Color(color[i][j]));
        return new Picture(picture);
    }

    // width of current picture
    public int width() {
        return width;
    }

    // height of current picture
    public int height() {
        return height;
    }

    // energy of pixel at column x and row y
    public double energy(int x, int y) {
        if (x < 0 || y < 0 || x >= width() || y >= height())
            throw new IndexOutOfBoundsException();

        return energy[y][x];
    }

    // sequence of indices for horizontal seam
    public int[] findHorizontalSeam() {

        distToSink = Double.POSITIVE_INFINITY;
        edgeToSink = Integer.MAX_VALUE;

        distTo = new double[height][width];
        edgeTo = new int[height][width];

        for (double[] r: distTo) Arrays.fill(r, Double.POSITIVE_INFINITY);
        for (int[] r: edgeTo) Arrays.fill(r, Integer.MAX_VALUE);

        for (int i = 0; i < height(); i++) {
            distTo[i][0] = 1000;
            edgeTo[i][0] = -1;
        }

        for (int depth = height() - 1; depth > 0; depth--)
            for (int out = 0; out < width() && depth + out < height(); out++)
                horizontalVisit(depth + out, out);

        for (int top = 0; top < width(); top++)
            for (int depth = 0; depth + top < width() && depth < height(); depth++)
                horizontalVisit(depth, depth + top);

        int[] seam = new int[width()];
        seam[width() - 1] = edgeToSink;
        for (int j = width() - 1; j > 0; j--)
            seam[j - 1] = edgeTo[seam[j]][j];

        distTo = null;
        edgeTo = null;

        return seam;
    }

    // sequence of indices for vertical seam
    public int[] findVerticalSeam() {

        distToSink = Double.POSITIVE_INFINITY;
        edgeToSink = Integer.MAX_VALUE;

        distTo = new double[height][width];
        edgeTo = new int[height][width];

        for (double[] r: distTo) Arrays.fill(r, Double.POSITIVE_INFINITY);
        for (int[] r: edgeTo) Arrays.fill(r, Integer.MAX_VALUE);

        Arrays.fill(distTo[0], 1000);
        Arrays.fill(edgeTo[0], -1);

        for (int top = width() - 1; top >= 0; top--)
            for (int depth = 0; depth < height() && top + depth < width(); depth++)
                verticalVisit(depth, depth + top);

        for (int depth = 1; depth < height(); depth++)
            for (int out = 0; depth + out < height() && out < width(); out++)
                verticalVisit(depth + out, out);

        int[] seam = new int[height()];
        seam[height() - 1] = edgeToSink;
        for (int i = height() - 1; i > 0; i--)
            seam[i - 1] = edgeTo[i][seam[i]];

        distTo = null;
        edgeTo = null;

        return seam;
    }

    private void horizontalVisit(int i, int j) {
        if (j == width() - 1) {
            relax(i, j, true);
        }
        else if (i == height() - 1) {
            relax(i, j, i, j + 1, true);
            relax(i, j, i - 1, j + 1, true);
        }
        else if (i == 0) {
            relax(i, j, i, j + 1, true);
            relax(i, j, i + 1, j + 1, true);
        }
        else {
            relax(i, j, i - 1, j + 1, true);
            relax(i, j, i, j + 1, true);
            relax(i, j, i + 1, j + 1, true);
        }
    }

    private void verticalVisit(int i, int j) {
        if (i == height() - 1) {
            relax(i, j, false);
        }
        else if (j == width() - 1) {
            relax(i, j, i + 1, j - 1, false);
            relax(i, j, i + 1, j, false);
        }
        else if (j == 0) {
            relax(i, j, i + 1, j, false);
            relax(i, j, i + 1, j + 1, false);
        }
        else {
            relax(i, j, i + 1, j - 1, false);
            relax(i, j, i + 1, j, false);
            relax(i, j, i + 1, j + 1, false);
        }
    }

    /**
     *
     * @param i - height
     * @param j - width
     * @param flag - show what type of removal we are doing: true - Horizontal, false - Vertical
     */
    private void relax(int i, int j, boolean flag) {
        if (validIndex(i, j)) {
            if (distToSink > distTo[i][j]) {
                distToSink = distTo[i][j];
                if (flag) edgeToSink = i;
                else edgeToSink = j;
            }
        }
    }

    /**
     *
     * @param i1 - height first pixel
     * @param j1 - width first pixel
     * @param i2 - height second pixel
     * @param j2 - width second pixel
     * @param flag - show what type of removal we are doing: true - Horizontal, false - Vertical
     */
    private void relax(int i1, int j1, int i2, int j2, boolean flag) {
        if (validIndex(i1, j1) && validIndex(i2, j2)) {
            if (distTo[i2][j2] > distTo[i1][j1] + energy[i2][j2]) {
                distTo[i2][j2] = distTo[i1][j1] + energy[i2][j2];
                if (flag) edgeTo[i2][j2] = i1;
                else edgeTo[i2][j2] = j1;
            }
        }
    }

    // Check that the index lies in the image itself
    private boolean validIndex(int i, int j) {
        return !(i < 0 || j < 0 || i > height() || j > width());
    }

    /**
     * Remove horizontal seam
     *
     * @param seam - the seam to be removed from the picture
     */
    public void removeHorizontalSeam(int[] seam) {
        if (height() <= 1)
            throw new IllegalArgumentException("Picture's height <= 1");
        if (seam == null)
            throw new NullPointerException();
        if (seam.length != width())
            throw new IllegalArgumentException("Argument is not correct");

        int yLast = seam[0];
        for (int y: seam) {
            if (y >= height() || y < 0)
                throw new IllegalArgumentException("Index out of bounds");
            if (Math.abs(y - yLast) > 1)
                throw new IllegalArgumentException("Index not adjacent");
            yLast = y;
        }

        int[][] colorUpd = new int[height() - 1][width()];
        double[][] energyUpd = new double[height() - 1][width()];

        for (int i = 0; i < width(); ++i) {

            for (int j = 0; j < seam[i]; ++j) colorUpd[j][i] = color[j][i];
            for (int j = seam[i] + 1; j < height(); ++j) colorUpd[j - 1][i] = color[j][i];
        }

        this.color = colorUpd;
        this.height--;

        for (int i = 0; i < height(); ++i) {
            for (int j = 0; j < width(); ++j) {
                energyUpd[i][j] = calculateEnergy(j, i);
            }
        }
        this.energy = energyUpd;

    }

    /**
     * Remove vertical seam
     *
     * @param seam - the seam to be removed from the picture
     */
    public void removeVerticalSeam(int[] seam) {
        if (width() <= 1)
            throw new IllegalArgumentException("Picture's height <= 1");
        if (seam == null)
            throw new NullPointerException();
        if (seam.length != height())
            throw new IllegalArgumentException("Argument is not correct");

        int yLast = seam[0];
        for (int y : seam) {
            if (y >= width() || y < 0)
                throw new IllegalArgumentException("Index out of bounds");
            if (Math.abs(y - yLast) > 1)
                throw new IllegalArgumentException("Index not adjacent");
            yLast = y;
        }

        int[][] colorUpd = new int[height()][width() - 1];
        double[][] energyUpd = new double[height()][width() - 1];

        for (int j = 0; j < height(); ++j) {

            for (int i = 0; i < seam[j]; ++i) colorUpd[j][i] = color[j][i];
            for (int i = seam[j] + 1; i < width(); ++i) colorUpd[j][i - 1] = color[j][i];
        }

        this.color = colorUpd;
        this.width--;

        for (int i = 0; i < height(); ++i) {
            for (int j = 0; j < width(); ++j) {
                energyUpd[i][j] = calculateEnergy(j, i);
            }
        }
        this.energy = energyUpd;

    }


    //  unit testing (optional)
    public static void main(String[] args) {
        
    }

}
