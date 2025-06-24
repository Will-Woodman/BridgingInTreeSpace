/**
    TreeViewer

    Copyright (C) 2011  Tom M. W. Nye

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact the author at:  <tom.nye@ncl.ac.uk>
                            <http://www.mas.ncl.ac.uk/~ntmwn/>
 */

package treegraphics;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import javax.swing.BorderFactory;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenuItem;
import javax.swing.JPanel;
import javax.swing.JPopupMenu;
import javax.swing.JScrollPane;
import javax.swing.JSlider;
import javax.swing.RepaintManager;
import javax.swing.SwingUtilities;
import javax.swing.ToolTipManager;
import javax.swing.WindowConstants;

/**
 * Window in which trees which implement DrawTree can be drawn.
 * Based on old "DendrogramViewer" class
 *
 * Instance variables:
 * The principal object is the tree being viewed.
 *
 * Methods:
 * - Print: allows the graphics to be printed
 * - refresh view: called when the tree has changed in some way
 * - set viewer parameters: called when the viewer configuration form is closed to update the scaling
 */



public class TreeViewer extends JPanel implements MouseListener, ActionListener {

    // Instance variables
    // ------------------

    protected DrawTree theTreePlotter;
    private Dimension area; //indicates area taken up by graphics
    protected DrawingPane drawingPane;
    private JPopupMenu theBackgroundMenu; // Items for the view menu
    private JMenuItem scaleItem;

    // Constants controlling appearance of the tree
    public final int BORDER_WIDTH = 10;

    // Variables controlling the time scale
    // ------------------------------------

    // Variables associated with the sliders
    static final int scaleMin = 0;  // value corresponding to left hand end of scale
    static final int scaleMax = 99; // value corresponding to right hand end of scale
    static final int tickSpace = 10;
    // These aren't stored in the dialog as they're needed in the viewer
    // in order to calculate scales


    // Constants relating to the linear scaling
    protected final int smallestTreeSize = 100;  // Smallest allowable tree size in pix. 250 is good.
    protected final int largestTreeSize = 2620; // Greatest size tree in pix
    // NB: this size has been determined to correspond to the edge of an A4 page
    // when the entire dendrogram is printed on one page.
    protected double pixPerUnitLength; // Converts time / branch length to pixels
    protected int scaleSliderValue = 30; // The initial value for the scale slider

        /** Creates a new instance of DendrogramViewer */
    public TreeViewer(DrawTree t) {
        super(new BorderLayout());
        theTreePlotter = t;

        // Set up the viewer parameters
        setViewerParameters(scaleSliderValue);

        // Prepare the drawing panel
        area = new Dimension();
        drawingPane = new DrawingPane();
        drawingPane.setBackground(Color.white);
        drawingPane.addMouseListener(this);
        ToolTipManager.sharedInstance().registerComponent(drawingPane);

        //Put the drawing area in a scroll pane.
        JScrollPane scroller = new JScrollPane();
        scroller.setViewportView(drawingPane);
        //Lay out
        add(scroller, BorderLayout.CENTER);

        // Set up the viewer menu
        theBackgroundMenu = new JPopupMenu();
        scaleItem = new JMenuItem("Configure scaling...");
        scaleItem.addActionListener(this);
        theBackgroundMenu.add(scaleItem);
        refreshView();
    }

    /** Refresh the drawing */
    public void refreshView() {
        if (theTreePlotter == null) {
            drawingPane.revalidate();
            drawingPane.repaint();
            return;
        }
        area = theTreePlotter.calculateSize(this, drawingPane.getGraphics());
        drawingPane.setPreferredSize(area);
        drawingPane.revalidate();
        drawingPane.repaint();
        drawingPane.requestFocusInWindow();
    }


    /** Change the tree eg when the user has loaded up a new one */
    public void changeTree(DrawTree t) {
        theTreePlotter = t;
        setViewerParameters(scaleSliderValue);
        refreshView();
    }

    /** Print -- see the drawingPane print method */
    public void print() {
        java.awt.print.PrinterJob printJob = java.awt.print.PrinterJob.getPrinterJob();
        javax.print.PrintService s = printJob.getPrintService();
        printJob.setPrintable(drawingPane);
        if (printJob.printDialog())
          try {
            printJob.print();
          } catch(java.awt.print.PrinterException pe) {
            System.out.println("Error printing: " + pe);
        }
    }

    public JPanel getDrawingPane() {
        return drawingPane;
    }


    // ----------------------------------------------------------------------

    /* Stuff to do with scaling and coordinates */

    /** Set the viewer parameters given slider values.
     This is called from the Viewer form and constructor to convert values from
     the sliders into useful parameter values. */
    public void setViewerParameters() {
        setViewerParameters(scaleSliderValue);
    }
    public void setViewerParameters(int scaleValue) {
        if (theTreePlotter == null) return; // Don't do anything if the tree is null!

        /* The size of the dendrogram depends in an exponential way on the slider
         ie. adjustment is finer for smaller branch lengths
         This is done so that the minimum tree width is smallestTreeSize, and
         the max largestTreeSize */

        double logymax = Math.log(largestTreeSize);
        double logymin = Math.log(smallestTreeSize);
        double a = (logymax-logymin)/(scaleMax-scaleMin);
        double x0 = (scaleMin*logymax - scaleMax*logymin)/(logymax-logymin);

        double s = Math.exp(a*(scaleValue-x0));
        pixPerUnitLength = s/theTreePlotter.getMaxPathLength();

    }
    public int getScaleValue() {
        return scaleSliderValue;
    }

    /* MATH METHODS for coordinates */

    /** Convert a branch length to a distance in pixels */
    public int convertBranchLenToPix(double dt) {
        int i;
        double d;
        d = dt * pixPerUnitLength;
        i = (int) Math.floor(d+0.5);
        return i;
    }

    public int polarToX(double direction, double dt) {
        int i;
        double d;
        d = dt * pixPerUnitLength * Math.cos(direction);
        i = (int) Math.floor(d+0.5);
        return i;
    }

    public int polarToY(double direction, double dt) {
        int i;
        double d;
        d = dt * pixPerUnitLength * Math.sin(direction);
        i = (int) Math.floor(d+0.5);
        return i;
    }

    // -------------------------------------------------------------------

    /* Deal with menus and mice */

    /** Cheeky menu stuff -- add an item.
     Useful for customizing the viewer. */
    public void addToBackgroundMenu(JMenuItem theItem) {
        theBackgroundMenu.add(theItem);
    }

    /** Handle events in the treeViewer -- response to background menu */
    public void actionPerformed(ActionEvent e) {
        JMenuItem source = (JMenuItem)(e.getSource());
        // Show the form if user clicked on the background menu
        if (source == scaleItem) {
            showScaleForm();
            drawingPane.repaint();
        }
    }

    /** Show the diaglog box to edit the view parameters */
    public void showScaleForm() {
         ScaleForm theForm = null;
         Container theContainer = this.getTopLevelAncestor();
         if (theContainer instanceof JFrame)
              theForm = new ScaleForm((JFrame) theContainer);
         else if (theContainer instanceof JDialog)
             theForm = new ScaleForm((JDialog) theContainer);
         else
             theForm = new ScaleForm(new JFrame());
         theForm.show();
    }


    /** Handle mouse down
     Right button down brings up popup menu on the tree */
    public void mousePressed(MouseEvent evt) {
        // Search to see whether the mouse was clicked on a node
        if (theTreePlotter != null) {
            treebase.Graph.Vertex v = theTreePlotter.getVertexAtPoint(evt.getPoint());
            if (v != null) {
                // Do nothing
            }
            else if (SwingUtilities.isRightMouseButton(evt)) {
                // Display a menu allowing access to the view form
                theBackgroundMenu.show(drawingPane, evt.getPoint().x, evt.getPoint().y);
            }
        }
    }

    public void mouseReleased(MouseEvent e){}
    public void mouseClicked(MouseEvent e){}
    public void mouseEntered(MouseEvent e){}
    public void mouseExited(MouseEvent e){}

    // -----------------------------------------------------------------------
    // SCLAEFORM inner classs

    /** ScaleForm is a dialog box that enables the user to edit the
     parameters controlling the view of the tree. */
    public class ScaleForm extends JDialog implements ActionListener {

        // Components
        private JSlider scaleSlider;

        /** Creates a new instance of ViewerForm. */
        public ScaleForm(JFrame parent) {
            super(parent, "TreeViewer Set-up", true);
            buildForm();
        }
        public ScaleForm(JDialog parent) {
            super(parent, "TreeViewer Set-up", true);
            buildForm();
        }

        /** Build the components in the form */
        private void buildForm() {

            JButton OKButton;
            JPanel buttonPanel;
            JPanel scalePanel;

            // Set up the scale slider
            scaleSlider = new JSlider(JSlider.HORIZONTAL, scaleMin, scaleMax, scaleSliderValue);
            scaleSlider.setMajorTickSpacing(tickSpace);
            scaleSlider.setPaintTicks(true);
            java.util.Hashtable scaleLabelTable = new java.util.Hashtable();
            scaleLabelTable.put(new Integer(scaleMin), new JLabel("Smaller"));
            scaleLabelTable.put(new Integer(scaleMax), new JLabel("Larger"));
            scaleSlider.setLabelTable(scaleLabelTable);
            scaleSlider.setPaintLabels(true);
            scalePanel = new JPanel();
            scalePanel.setLayout(new FlowLayout(FlowLayout.CENTER));
            scalePanel.setBorder(BorderFactory.createCompoundBorder(
            BorderFactory.createTitledBorder("Scale"),
            BorderFactory.createEmptyBorder(5,5,5,5)));
            scalePanel.add(scaleSlider);

            buttonPanel = new JPanel();
            buttonPanel.setLayout(new FlowLayout(FlowLayout.RIGHT));
            buttonPanel.setMaximumSize(new Dimension(32767, 200));
            OKButton = new JButton();
            OKButton.setText("OK");
            OKButton.setActionCommand("OK");
            OKButton.addActionListener(this);

            buttonPanel.add(OKButton);

            getContentPane().setLayout(new BoxLayout(getContentPane(), BoxLayout.Y_AXIS));
            setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
            getContentPane().add(scalePanel);
            getContentPane().add(buttonPanel);
            getRootPane().setDefaultButton(OKButton);

            pack();
        }

        /** Handle a click on ok or the radio buttons */
        public void actionPerformed(ActionEvent e) {
            // Handle a click on OK
            if ("OK".equals(e.getActionCommand())) {
                // Recalculate view parameters
                scaleSliderValue = scaleSlider.getValue();
                TreeViewer.this.setViewerParameters(scaleSliderValue);
                this.dispose();
                // Refresh the treeViewer
                TreeViewer.this.refreshView();
            }
        }
        // End the ScaleForm inner class
    }

    // End form showing scaler
    // ----------------------------------------------------------------------


    // ----------------------------------------------------------------------
    // Start DrawingPane inner class

    /** DrawingPane is the graphics pane that lies inside the TreeViewers's scroll pane
     It handles the drawing of the tree */
    protected class DrawingPane extends JPanel implements java.awt.print.Printable {

        // Override the default paint method: draw the tree!
        protected void paintComponent(Graphics g) {
            super.paintComponent(g);
            if (theTreePlotter == null) return; // Don't do anything if the tree is null!
            Dimension treeArea = new Dimension();
            theTreePlotter.draw(g, TreeViewer.this, treeArea);
            area.width = treeArea.width;
            area.height = treeArea.height;
            this.setPreferredSize(area);
        }

        // Display tool tip information: this is associated to branches rather than nodes
        public String getToolTipText(MouseEvent evt) {
            if (theTreePlotter != null) {
                treebase.Graph.Edge e = theTreePlotter.getEdgeAtPoint(evt.getPoint());
                if (e != null) {
                    String s = theTreePlotter.getEdgeLabel(e);
                    if (s.isEmpty()) return null;
                    else return s;
                }
                else {
                    treebase.Graph.Vertex v = theTreePlotter.getVertexAtPoint(evt.getPoint());
                    if (v != null) {
                        String s = theTreePlotter.getVertexLabel(v);
                        if (s.isEmpty()) return null;
                        else return s;
                    }
                }
            }
            return null;
        }

        /** Print the pane
         The dendrogram is scaled to fill a single page of A4 */
        public int print(Graphics g, java.awt.print.PageFormat pageFormat, int pageIndex) throws java.awt.print.PrinterException {
            if (pageIndex > 0) {
                return(NO_SUCH_PAGE);
            }
            else {
                Graphics2D g2d = (Graphics2D)g;
                g2d.translate(pageFormat.getImageableX(), pageFormat.getImageableY());
                // Turn off double buffering
                RepaintManager currentManager = RepaintManager.currentManager(this);
                currentManager.setDoubleBufferingEnabled(false);

                Dimension d = theTreePlotter.calculateSize(TreeViewer.this, g2d);
                double fac;
                double r = ((double) d.height)/((double) d.width);
                if (r>1.41) {
                    // image is too tall for page
	            fac = 9.50/(d.height/72.0);
                }
                else {
                    // image is to wide for page
                    fac = 430.0/d.width;
                }
                g2d.scale(fac,fac);
                Dimension treeArea = new Dimension();
                theTreePlotter.draw(g2d, TreeViewer.this, treeArea);

                // Turn double buffering back on
                currentManager.setDoubleBufferingEnabled(true);
                return(PAGE_EXISTS);
            }
        }

    }

    // End DrawingPane inner class
    // ------------------------------------------------------------------------


    /* Testing area */

    public static void main(String[] args) throws treebase.AlgorithmError {

        treebase.RootedTree t = null;
        JFileChooser fc = new JFileChooser();
        int returnVal = fc.showOpenDialog(null);
        if (returnVal == JFileChooser.APPROVE_OPTION) {
            java.io.File theFile = fc.getSelectedFile();
            try {
                t = new treebase.RootedTree(theFile);
            }
            catch (Exception problemMakingTestTree) {
                System.out.println("IOException. "+problemMakingTestTree.getMessage());
            }
        }

        TreePlotter thePlotter = new TreePlotter(t);

        //Create and set up window.
        JFrame frame = new JFrame("TreeViewer Demo");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

        // Create the treeViewer
        TreeViewer theViewer = new TreeViewer(thePlotter);

        //Set up the content pane.
        theViewer.setOpaque(true); //content panes must be opaque
        frame.setContentPane(theViewer);

        // Set position and size
        frame.setLocation(100,100);

        //Display the window.
        frame.pack();
        theViewer.refreshView();
        frame.setVisible(true);

        // Print
//        theViewer.print();

    }

}
