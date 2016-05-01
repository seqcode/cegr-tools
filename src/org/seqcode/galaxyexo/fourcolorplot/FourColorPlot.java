package org.seqcode.galaxyexo.fourcolorplot;

/**
 * FourColorPlot: 
 * 
 * Produces a graphical representation of FASTA data with each nucleotide represented by a selected color.
 * 
 * @author Will Lai
 */
import java.awt.Color;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.sql.Timestamp;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;
import java.util.Scanner;

import javax.imageio.ImageIO;

public class FourColorPlot {
	private static File INPUT = null;
	private static File OUTPUT = null;
	
	private static Color AColor = new Color(254, 25, 24);
	private static Color TColor = new Color(50, 204, 60);
	private static Color GColor = new Color(252, 252, 80);
	private static Color CColor = new Color(43, 49, 246);
	private static Color NColor = Color.GRAY;
	
	private static int PIXEL_HEIGHT = 2;
	private static int PIXEL_WIDTH = 2;
	
	public static void main(String[] args) throws IOException {
		loadConfig(args);
		System.out.println("\n" + getTimeStamp());
	
		generatePLOT();
		
		System.out.println(getTimeStamp());
	}
	
	public static void generatePLOT() throws IOException {
		List<String> seq = new ArrayList<String>();
		int maxLen = 0;

		Scanner scan = new Scanner(INPUT);
		while (scan.hasNextLine()) {
			String temp = scan.nextLine();
			if(!temp.contains(">")) {
				if (maxLen < temp.length()) maxLen = temp.length();
				seq.add(temp);
			}
		}
		scan.close();
		int pixwidth = maxLen * PIXEL_WIDTH;
		int pixheight = seq.size() * PIXEL_HEIGHT;
		
		System.setProperty("java.awt.headless", "true");
		BufferedImage im = new BufferedImage(pixwidth, pixheight, BufferedImage.TYPE_INT_ARGB);
        Graphics g = im.getGraphics();
        Graphics2D g2 = (Graphics2D)g;
        g2.setColor(Color.WHITE);
        g2.fillRect(0, 0, pixwidth, pixheight);
        
        int count = 0;
        for (int x = 0; x < seq.size(); x++){
        	String s = seq.get(x);
        	char[] letters = s.toCharArray();
        	for (int y = 0; y < letters.length; y++){
        		switch(letters[y]){
        		case 'A':
        		case 'a':
        			g.setColor(AColor);
        			break;
        		case 'T':
        		case 't':
                    g.setColor(TColor);
        			break;
        		case 'G':
        		case 'g':
                    g.setColor(GColor);
        			break;
        		case 'C':
        		case 'c':
                    g.setColor(CColor);
        			break;
        		case '-':
                    g.setColor(Color.WHITE);
        			break;
                default:
                	g.setColor(NColor);
        		}
                g.fillRect(y * PIXEL_WIDTH, count * PIXEL_HEIGHT, PIXEL_WIDTH, PIXEL_HEIGHT);
        	}
            count++;
        }
        
        try { ImageIO.write(im, "png", OUTPUT);
        } catch (IOException ex) { ex.printStackTrace(); }
	}
	
	public static void loadConfig(String[] command){
		String A = null, T = null, G = null, C = null, N = null;
		for (int i = 0; i < command.length; i++) {
			switch (command[i].charAt((1))) {
				case 'f':
					INPUT = new File(command[i + 1]);
					i++;
					break;
				case 'o':
					OUTPUT = new File(command[i + 1]);
					i++;
					break;
				case 'h':
					PIXEL_HEIGHT = Integer.parseInt(command[i + 1]);
					i++;
					break;
				case 'w':
					PIXEL_WIDTH = Integer.parseInt(command[i + 1]);
					i++;
					break;
        		case 'A':
        		case 'a':
					A = command[i + 1];
					i++;
					break;
        		case 'T':
        		case 't':
					T = command[i + 1];
					i++;
					break;
        		case 'G':
        		case 'g':
					G = command[i + 1];
					i++;
					break;
        		case 'C':
        		case 'c':
					C = command[i + 1];
					i++;
					break;
        		case 'N':
        		case 'n':
					N = command[i + 1];
					i++;
					break;
			}
		}
		if(INPUT == null) {
			System.out.println("No FASTA File loaded!!!\n");
			printUsage();
			System.exit(1);
		}
		if(PIXEL_HEIGHT < 0 || PIXEL_WIDTH < 0) {
			System.out.println("Invalid Pixel Size Selected!!!\n");
			printUsage();
			System.exit(1);	
		}
		
		if(A != null) {
			Color tempColor = verifyRGB(A);
			if(tempColor == null) {
				System.out.println("Invalid A RGB Color!!!\n");
				printUsage();
				System.exit(1);
			} else { AColor = tempColor; }
		}
		if(T != null) {
			Color tempColor = verifyRGB(T);
			if(tempColor == null) {
				System.out.println("Invalid T RGB Color!!!\n");
				printUsage();
				System.exit(1);
			} else { TColor = tempColor; }
		}
		if(G != null) {
			Color tempColor = verifyRGB(G);
			if(tempColor == null) {
				System.out.println("Invalid G RGB Color!!!\n");
				printUsage();
				System.exit(1);
			} else { GColor = tempColor; }
		}
		if(C != null) {
			Color tempColor = verifyRGB(C);
			if(tempColor == null) {
				System.out.println("Invalid C RGB Color!!!\n");
				printUsage();
				System.exit(1);
			} else { CColor = tempColor; }
		}
		if(N != null) {
			Color tempColor = verifyRGB(N);
			if(tempColor == null) {
				System.out.println("Invalid N RGB Color!!!\n");
				printUsage();
				System.exit(1);
			} else { NColor = tempColor; }
		}
				
		if(OUTPUT == null) {
			OUTPUT = new File(INPUT.getName().split("\\.")[0] + ".png");
		}
				
		System.out.println("-----------------------------------------\nCommand Line Arguments:");
		System.out.println("FASTA file: " + INPUT);
		System.out.println("Output: " + OUTPUT);
		System.out.println("Pixel Height per Nucleotide: " + PIXEL_HEIGHT);
		System.out.println("Pixel Width per Nucleotide: " + PIXEL_WIDTH);
		System.out.println("A Color: " + AColor.toString());
		System.out.println("T Color: " + TColor.toString());
		System.out.println("G Color: " + GColor.toString());
		System.out.println("C Color: " + CColor.toString());
		System.out.println("N Color: " + NColor.toString());
	}
	
	public static Color verifyRGB(String color) {
		Color newColor = null;
		if(color == null) { return newColor; }
		String[] rgb = color.split(",");
		if(rgb.length != 3) { return newColor; }
		int RED = Integer.parseInt(rgb[0]);
		if(RED < 0 || RED > 255) { return newColor; }
		int GREEN = Integer.parseInt(rgb[1]);
		if(GREEN < 0 || GREEN > 255) { return newColor; }
		int BLUE = Integer.parseInt(rgb[2]);
		if(BLUE < 0 || BLUE > 255) { return newColor; }
		newColor = new Color(RED, GREEN, BLUE);
		return newColor;
	}
		
	public static void printUsage() {
		System.err.println("Usage: java -jar FourColorPlot.jar -f [FASTAFile] [Options]");
		System.err.println("-----------------------------------------");
		System.err.println("\nRequired Parameter:");
		System.err.println("FASTA File:\t\t\t-f\tFASTA file");
		
		System.err.println("\nOptional Parameters:");
		System.err.println("Output File Name:\t\t-o\tOutput file");
		System.err.println("Pixel Height per Nucleotide:\t-h\t(Default = 2)");
		System.err.println("Pixel Width per Nucleotide:\t-w\t(Default = 2)");
		System.err.println("A Color RGB comma-delimited:\t-A\t(Default 254,25,24)");
		System.err.println("T Color RGB comma-delimited:\t-T\t(Default 50,204,60)");
		System.err.println("G Color RGB comma-delimited:\t-G\t(Default 252,252,80)");
		System.err.println("C Color RGB comma-delimited:\t-C\t(Default 43,49,246)");
		System.err.println("N Color RGB comma-delimited:\t-N\t(Default 128,128,128)");
	}
	
	private static String getTimeStamp() {
		Date date= new Date();
		String time = new Timestamp(date.getTime()).toString();
		return time;
	}
}
