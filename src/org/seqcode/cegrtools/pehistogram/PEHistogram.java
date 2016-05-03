package org.seqcode.cegrtools.pehistogram;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;

import htsjdk.samtools.AbstractBAMFileIndex;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloseableIterator;

/**
 * PEHistogram
 * 
 * Creates a histogram of inferred fragment size lengths from mapped paired-end sequencing data. 
 * 
 * @author Will Lai
 *
 */
public class PEHistogram {
	public File PATH = null;
	
	public static File BAM = null;
	public static File BAI = null;
	public static String OUTPUT = null;
	
	public static int MinSize = 0;
	public static int MaxSize = 1000;
			
	public static void main(String[] args) {
		loadConfig(args);
		System.out.println("Generating Statistics...");
		generateHist();
		System.out.println("Program Complete");
	}
	
	public static void generateHist() {
		PrintStream OUT = null;
		try {
			OUT = new PrintStream(OUTPUT + ".out");
		} catch (FileNotFoundException e1) {
			e1.printStackTrace();
		}
		if(OUT != null) OUT.println(BAM.getName());
		if(OUT != null) OUT.println("Chromosome_ID\tChromosome_Size\tAligned_Reads");
		
		//Code to get individual chromosome stats
		SAMFileReader reader = new SAMFileReader(BAM, BAI);
		AbstractBAMFileIndex bai = (AbstractBAMFileIndex) reader.getIndex();
		double totalTags = 0;
		double totalGenome = 0;
	
		for (int z = 0; z < bai.getNumberOfReferences(); z++) {
			SAMSequenceRecord seq = reader.getFileHeader().getSequence(z);
			double aligned = reader.getIndex().getMetaData(z).getAlignedRecordCount();
			if(OUT != null) OUT.println(seq.getSequenceName() + "\t" + seq.getSequenceLength() + "\t" + aligned);
			totalTags += aligned;
			totalGenome += seq.getSequenceLength();
		}
		if(OUT != null) OUT.println("Total Genome Size: " + totalGenome + "\tTotal Aligned Tags: " + totalTags + "\n");
		
		//Output replicates used to make bam file
		for( String comment : reader.getFileHeader().getComments()) {
			if(OUT != null) OUT.println(comment);
		}
		
		//Output program used to align bam file
		for (int z = 0; z < reader.getFileHeader().getProgramRecords().size(); z++) {
			if(OUT != null) {
				OUT.print(reader.getFileHeader().getProgramRecords().get(z).getId() + "\t");
				OUT.println(reader.getFileHeader().getProgramRecords().get(z).getProgramVersion());
				OUT.println(reader.getFileHeader().getProgramRecords().get(z).getCommandLine());
			}
		}		
		if(OUT != null) OUT.println();
		
		double average = 0;
		double counter = 0;
		double[] HIST = new double[(MaxSize - MinSize) + 1];
		
		CloseableIterator<SAMRecord> iter = reader.iterator();
		while (iter.hasNext()) {
			SAMRecord sr = iter.next();
			if(sr.getReadPairedFlag()) {
				if(sr.getProperPairFlag() && sr.getFirstOfPairFlag()) {
					int distance = Math.abs(sr.getInferredInsertSize());
					if(distance <= MaxSize && distance >= MinSize) HIST[distance - MinSize]++;
					average += distance;
					counter++;
				}
			}
		}
		iter.close();
		if(counter != 0) average /= counter;
		
		if(OUT != null) OUT.println("Average Insert Size: " + average + "\nNumber of ReadPairs: " + counter + "\n\nHistogram\nSize (bp)\tFrequency");
		int[] DOMAIN = new int[(MaxSize - MinSize) + 1];
		for(int z = 0; z < HIST.length; z++) {
			int bp = MinSize + z;
			DOMAIN[z] = bp;
			if(OUT != null) OUT.println(bp + "\t" + HIST[z]);
		}
		reader.close();
		bai.close();
		
		try {
			Histogram.createBarChart(OUTPUT, HIST, DOMAIN);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public static void loadConfig(String[] command){
		for (int i = 0; i < command.length; i++) {
			switch (command[i].charAt((1))) {
				case 'B':
					BAM = new File(command[i + 1]);
					i++;
					break;
				case 'l':
					MinSize = Integer.parseInt(command[i + 1]);
					i++;
					break;
				case 'u':
					MaxSize = Integer.parseInt(command[i + 1]);
					i++;
					break;
				case 'o':
					OUTPUT = command[i + 1];
					i++;
					break;
			}
		}
		if(BAM == null) {
			System.err.println("No BAM File Loaded!!!\n");
			printUsage();
			System.exit(0);
		}
		BAI = new File(BAM.getAbsolutePath() + ".bai");
		if(!BAI.exists() || BAI.isDirectory()) {
			System.err.println("No BAI File Detected!!!\n");
			printUsage();
			System.exit(1);
		}
		if(OUTPUT == null) {
			OUTPUT = BAM.getName().split("\\.")[0] + "_STATS";
		}
		System.out.println("-----------------------------------------\nCommand Line Arguments:");
		System.out.println("BAM file: " + BAM.getName());
		System.out.println("BAI file: " + BAI.getName());
		System.out.println("Output filename: " + OUTPUT);
		System.out.println("-----------------------------------------\n");

	}

	public static void printUsage() {
		System.err.println("Usage: java -jar PEHistogram.jar -B [BAMFile] [options]");
		System.err.println("-----------------------------------------");
		System.err.println("Required Parameter:");
		System.err.println("BAM File: -B [BAM File]");
		System.err.println("Optional Parameters:");
		System.err.println("Set Lower Limit: -l [Lower Limit]");
		System.err.println("Set Upper Limit: -u [Upper Limit]");
		System.err.println("Set Output Filename:  -o [OutputFile]");
	}
}
