package org.seqcode.galaxyexo.pairedendcrossplot;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.sql.Timestamp;
import java.util.ArrayList;
import java.util.Date;
import java.util.Scanner;

/**
 * PairedendCrossPlot
 * 
 * Creates a cross-plot from paired end sequencing data. 
 * 
 * @author Will Lai
 *
 */
public class PairedendCrossPlot {
	
	private static File READ1 = null;
	private static File READ2 = null;
	private static File OUTPUT = null;
	private static int WINDOW = 500;
	
	private static int[] BP = null;
	private static int[] FOR = null;
	private static int[] REV = null;
	
	private static ArrayList<String> chrName = null;
	private static ArrayList<Integer> chrStart = null;
	
	private static int[] FHIST = null;
	private static int[] RHIST = null;
	
	public static void main(String[] args) throws FileNotFoundException {
		loadConfig(args); //Load configuration data
		FHIST = new int[(WINDOW * 2) + 1]; //Initialize master histograms
		RHIST = new int[(WINDOW * 2) + 1]; //Initialize master histograms
		System.out.println("\n" + getTimeStamp());
	
		//Pass through the Read2 file once in order to identify the start and stops of all chromosomes present in file
		System.out.println("Indexing Read 2 File...");
		indexFile(READ2);
		System.out.println("Indexing Complete\n" + getTimeStamp());
	
		//Parse read 1
		System.out.println("\n" + getTimeStamp() + "\nParsing Read 1 File...");
		parseRead1(READ1);
		System.out.println("\nRead 1 Parsed\n" + getTimeStamp());
		
		//Output histogram
		System.out.println("\nOutputing final histogram...");
		outputHist(OUTPUT);
		System.out.println("Program Complete\n" + getTimeStamp());
		
	}
	
	public static void outputHist(File out) throws FileNotFoundException {
		PrintStream OUT = new PrintStream(out);
		for(int x = WINDOW * -1; x <= WINDOW; x++) {
			OUT.print("\t" + x);			
		}
		OUT.print("\nForward_Tags");
		for(int x = 0; x < FHIST.length; x++) {
			OUT.print("\t" + FHIST[x]);
		}
		OUT.print("\nReverse_Tags");
		for(int x = 0; x < RHIST.length; x++) {
			OUT.print("\t" + RHIST[x]);
		}
		OUT.close();
	}
	
	public static void parseRead1(File input) throws FileNotFoundException {
		Scanner scan = new Scanner(READ1);
		String currentChrom = "";
		int currentIndex = 0;
		while (scan.hasNextLine()) {
			String temp = scan.nextLine();
			if(!temp.contains("index") && !temp.contains("#")) {
				String[] array = temp.split("\t");
				int POS = Integer.parseInt(array[1]);
				int F = Integer.parseInt(array[2]);
				int R = Integer.parseInt(array[3]);
							
				//Check if newline's chrom is equal to current chrom in memory, reload new chrom if unequal
				if(!currentChrom.equals(array[0])) {
					BP = null; //reset BP immediately
					int newIndex = chrName.indexOf(array[0]);
					if(newIndex >= 0) {
						System.out.println("\nLoading: " + chrName.get(newIndex));
						loadRead2(READ2, chrName.get(newIndex), chrStart.get(newIndex), chrStart.get(newIndex + 1));
						currentChrom = array[0];
						currentIndex = 0;
						System.out.println("Loading Complete\n" + getTimeStamp());
					}
				}
				if(BP != null) {
					//Update local indexes to only upload data from local tag region defined by window
					while(POS - WINDOW > BP[currentIndex] && currentIndex < BP.length - 1) { currentIndex++; } //Find furthest upstream index in read2
					int currentStop = currentIndex;
			        while(POS + WINDOW > BP[currentStop] && currentStop < BP.length - 1) { currentStop++; } //Find furthest downstream index in read2
			        
			        //Populate histogram based on read2 tags flanking current position
			        int[] temp_F = new int[(WINDOW * 2) + 1];
			        int[] temp_R = new int[(WINDOW * 2) + 1];
			        for(int x = currentIndex; x <= currentStop; x++) {
			        	int index = BP[x] - POS + WINDOW;
			        	if(index >= 0 && index < temp_F.length) {
			        		temp_F[index] += FOR[x];
			        		temp_R[index] += REV[x];
			        	}
			        }
					
			        //Populate master histogram based on score
			        for(int x = 0; x < temp_F.length; x++) {
			                FHIST[x] += (temp_F[x] * F);
			                RHIST[x] += (temp_R[x] * F);
			        }
			        //Reverse arrays
			        for(int x = 0; x < FHIST.length / 2; x++) {
			        	int tempF = temp_F[x];
			        	int tempR = temp_R[x];
			        	temp_F[x] = temp_F[temp_F.length - x - 1];
			            temp_F[temp_F.length - x - 1] = tempF;
			            temp_R[x] = temp_R[temp_R.length - x - 1];
			            temp_R[temp_R.length - x - 1] = tempR;
			        }
			        for(int x = 0; x < FHIST.length; x++) {
			                FHIST[x] += (temp_R[x] * R);
			                RHIST[x] += (temp_F[x] * R);
			        }
				}
			}	
		}
		scan.close();
	}
	
	private static void loadRead2(File input, String currentChrom, int lineStart, int lineStop) throws FileNotFoundException {
		BP = new int[lineStop - lineStart - 1];
		FOR = new int[lineStop - lineStart - 1];
		REV = new int[lineStop - lineStart - 1];
		
		int counter = 0;
		Scanner scan = new Scanner(READ2);
		while (scan.hasNextLine()) {
			String temp = scan.nextLine();
				if(counter > lineStart) { 
					String[] array = temp.split("\t");
					if(array[0].equals(currentChrom)) {
						BP[counter - lineStart - 1] = Integer.parseInt(array[1]);
						FOR[counter - lineStart - 1] = Integer.parseInt(array[2]);
						REV[counter - lineStart - 1] = Integer.parseInt(array[3]);
					}
				}
			counter++;
		}
		scan.close();
		//System.out.println(BP.length + "\t" + FOR.length + "\t" + REV.length);
	}
	
	//Index Read 2 file for faster processing of chromosomes
	public static void indexFile(File input) throws FileNotFoundException {
		chrName = new ArrayList<String>();
		chrStart = new ArrayList<Integer>();
		Scanner scan = new Scanner(input);
		int linenumber = 0;
		while (scan.hasNextLine()) {
			String temp = scan.nextLine();
			if(!temp.contains("index") && !temp.contains("#")) {
				String[] array = temp.split("\t");
				if(!chrName.contains(array[0])) {
					chrName.add(array[0]);
					chrStart.add(new Integer(linenumber));
				}
			}
			linenumber++;
		}
		chrStart.add(new Integer(linenumber));
		scan.close();
	}
	
	public static void loadConfig(String[] command){
		for (int i = 0; i < command.length; i++) {
			switch (command[i].charAt((1))) {
				case 'f':
					READ1 = new File(command[i + 1]);
					i++;
					break;
				case 'r':
					READ2 = new File(command[i + 1]);
					i++;
					break;
				case 'o':
					OUTPUT = new File(command[i + 1]);
					i++;
					break;
				case 'w':
					WINDOW = Integer.parseInt(command[i + 1]);
					i++;
					break;
			}
		}
		if(READ1 == null) {
			System.out.println("No Read 1 File loaded!!!\n");
			printUsage();
			System.exit(0);
		}
		if(READ2 == null) {
			System.out.println("No Read 2 File loaded!!!\n");
			printUsage();
			System.exit(0);
		}
		if(WINDOW < 1) {
			System.out.println("Invalid Window Size Selected!!!\n");
			printUsage();
			System.exit(0);	
		}
						
		if(OUTPUT == null) {
			OUTPUT = new File(READ1.getName().split("\\.")[0] + ".out");
		}
				
		System.out.println("-----------------------------------------\nCommand Line Arguments:");
		System.out.println("Read 1 file: " + READ1);
		System.out.println("Read 2 file: " + READ2);
		System.out.println("Output: " + OUTPUT);
		System.out.println("Window Size: " + WINDOW);
	}
		
	public static void printUsage() {
		System.out.println("Usage: java -jar PairedendCrossPlot.jar -f [Read1File] -r [Read2File] [Options]");
		System.out.println("-----------------------------------------");
		System.out.println("\nRequired Parameter:");
		System.out.println("Read 1 File:\t\t\t-f\tRead 1 scIDX file");
		System.out.println("Read 2 File:\t\t\t-r\tRead 1 scIDX file");
		System.out.println("\nOptional Parameters:");
		System.out.println("Output File Name:\t\t-o\tOutput file");
		System.out.println("Window Size +/- Reference:\t-w\t(Default = 500)");
	}
	
	private static String getTimeStamp() {
		Date date= new Date();
		String time = new Timestamp(date.getTime()).toString();
		return time;
	}
	
}
