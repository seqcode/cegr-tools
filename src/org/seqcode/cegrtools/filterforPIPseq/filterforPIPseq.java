package org.seqcode.cegrtools.filterforPIPseq;

import htsjdk.samtools.AbstractBAMFileIndex;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.sql.Timestamp;
import java.util.ArrayList;
import java.util.Date;

public class filterforPIPseq {
	private static File bamFile = null;
	private static File genome = null;
	private static File output = null;
	private static String SEQ = "";
	
	static boolean FASTA_INDEX = true;
		
	public static void main(String[] args) throws IOException, InterruptedException {
		loadConfig(args); //Load configuration
		
		System.out.println("\n" + getTimeStamp());
		
		processREADS(); //Begin processing reads in BAM file
		
		System.out.println(getTimeStamp());
	}
	
	public static void processREADS() throws IOException, InterruptedException {
		File FAI = new File(genome + ".fai");
		//Check if FAI index file exists
		if(!FAI.exists() || FAI.isDirectory()) {
			System.out.println("FASTA Index file not found.\nGenerating new one...\n");
			FASTA_INDEX = buildFASTAIndex(genome);
		}		
		
		//Check if BAI index file exists
		File f = new File(bamFile + ".bai");
		if(f.exists() && !f.isDirectory()) {
			IndexedFastaSequenceFile QUERY = new IndexedFastaSequenceFile(genome);
			
			IOUtil.assertFileIsReadable(bamFile);
			IOUtil.assertFileIsWritable(output);
			final SamReader reader = SamReaderFactory.makeDefault().open(bamFile);
			reader.getFileHeader().setSortOrder(SAMFileHeader.SortOrder.coordinate);
			final SAMFileWriter writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(reader.getFileHeader(), false, output);
						
			//Code to get individual chromosome stats
			AbstractBAMFileIndex bai = (AbstractBAMFileIndex) reader.indexing().getIndex();
			for (int z = 0; z < bai.getNumberOfReferences(); z++) {
				SAMSequenceRecord seq = reader.getFileHeader().getSequence(z);
				System.out.println(seq.getSequenceName());
						
				CloseableIterator<SAMRecord> iter = reader.query(seq.getSequenceName(), 0, seq.getSequenceLength(), false);
				while (iter.hasNext()) {
					//Create the record object 
					SAMRecord sr = iter.next();
					if(sr.getReadPairedFlag()) {
						if(sr.getProperPairFlag() && sr.getFirstOfPairFlag()) {
							String filter = "";
							//if on the positive strand
							if(!sr.getReadNegativeStrandFlag()) {
								if(sr.getUnclippedStart() - 1 > 0) { filter = new String(QUERY.getSubsequenceAt(seq.getSequenceName(), sr.getUnclippedStart() - 1, sr.getUnclippedStart() - 1).getBases()); }
							}
							else {
								if(sr.getUnclippedEnd() + 1 <= seq.getSequenceLength()) {
										filter = new String(QUERY.getSubsequenceAt(seq.getSequenceName(), sr.getUnclippedEnd() + 1, sr.getUnclippedEnd() + 1).getBases());
										filter = RevComplement(filter);
								}
							}
							//System.out.println(sr.getReadString() + "\t" + seq.getSequenceName() + "\t" + sr.getUnclippedStart() + "\t" + sr.getUnclippedEnd() + "\t" + sr.getReadNegativeStrandFlag() + "\t" + filter);
							if(filter.toUpperCase().equals(SEQ)) { writer.addAlignment(sr); }							
						}
					} else {
						String filter = "";
						//if on the positive strand
						if(!sr.getReadNegativeStrandFlag()) {
							filter = new String(QUERY.getSubsequenceAt(seq.getSequenceName(), sr.getUnclippedStart() - 1, sr.getUnclippedStart() - 1).getBases());
						}
						else {
							filter = new String(QUERY.getSubsequenceAt(seq.getSequenceName(), sr.getUnclippedEnd() + 1, sr.getUnclippedEnd() + 1).getBases());
							filter = RevComplement(filter);
						}
						//System.out.println(sr.getReadString() + "\t" + seq.getSequenceName() + "\t" + sr.getUnclippedStart() + "\t" + sr.getUnclippedEnd() + "\t" + sr.getReadNegativeStrandFlag() + "\t" + filter);
						if(filter.toUpperCase().equals(SEQ)) { writer.addAlignment(sr); }		
					}
				}
				iter.close();
			}
			QUERY.close();
			writer.close();
			reader.close();
			bai.close();
		} else {
			System.out.println("BAI Index File does not exist for: " + bamFile.getName());
		}
	}

	/*
	 * Adapted from:
	 * https://github.com/mdshw5/pyfaidx/blob/master/pyfaidx/__init__.py
	 * pyfaidx python program for manipulating fasta files efficiently
	 *
	 *	contig_name\tcontig_length\toffset_distance_from_last_contig\tcolumnlength\tcolumnlength_with_endline\n"
     *	chr1    230218  6       60      61 
     *	chr2    813184  234067  60      61 
     */
    public static boolean buildFASTAIndex(File fasta) throws IOException {  	
    	boolean properFASTA = true;
    	ArrayList<String> IMPROPER_FASTA = new ArrayList<String>();
    	int counter = 0;

    	String contig = "";
    	int binaryOffset = 0;
    	int currentOffset = 0;
    	int contigLength = 0;
    	int column_Length = 0;
    	int untrimmed_Column_Length = 0;
    	    	
    	BufferedReader b_read = new BufferedReader(new FileReader(fasta));
    	//LineReader reader = new LineReader(b_read);
    	PrintStream FAI = new PrintStream(fasta.getCanonicalPath() + ".fai");
    	
    	String strLine = "";
    	while(!(strLine = readLine(b_read)).equals("")) {
    		//Pull parameters line
    		int current_untrimmed_Column_Length = strLine.length();
			int current_column_Length = strLine.trim().length();

			if(strLine.contains(">")) {
				if(IMPROPER_FASTA.size() > 1) {
					System.out.println("Unequal column size FASTA Line at:");
					for(int z = 0; z < IMPROPER_FASTA.size(); z++) {	System.out.println(contig + "\t" + IMPROPER_FASTA.get(z));	}
					properFASTA = false;
					break;
				}
				if(counter > 0) { FAI.println(contig + "\t" + contigLength + "\t" + currentOffset + "\t" + column_Length + "\t" + untrimmed_Column_Length);	}
				//Reset parameters for new contig
				untrimmed_Column_Length = 0;
				contigLength = 0;
				column_Length = 0;
				contig = strLine.trim().substring(1);
				binaryOffset += current_untrimmed_Column_Length;
				currentOffset = binaryOffset;
				IMPROPER_FASTA = new ArrayList<String>();
			} else {
				if(untrimmed_Column_Length == 0) { untrimmed_Column_Length = current_untrimmed_Column_Length; }
				if(column_Length == 0) { column_Length = current_column_Length;	}
				binaryOffset += current_untrimmed_Column_Length;
				contigLength += current_column_Length;
				
				//Check to make sure all the columns are equal. Index is invalid otherwise
				if(current_untrimmed_Column_Length != untrimmed_Column_Length || current_untrimmed_Column_Length == 0) { IMPROPER_FASTA.add(strLine.trim());	}
			}
			counter++;
    	}
		FAI.println(contig + "\t" + contigLength + "\t" + currentOffset + "\t" + column_Length + "\t" + untrimmed_Column_Length);
		b_read.close();
    	FAI.close();
    	
		if(properFASTA) System.out.println("Genome Index Built");
		else { new File(fasta.getName() + ".fai").delete(); }
		
		return properFASTA;
    }
    
	public static String RevComplement(String SEQ) {
		SEQ = SEQ.toUpperCase();
		String RC = "";
		for (int x = 0; x < SEQ.length(); x++){
			if(SEQ.charAt(x) == 'A') { RC = 'T' + RC; }
			else if(SEQ.charAt(x) == 'T') { RC = 'A' + RC; }
			else if(SEQ.charAt(x) == 'G') { RC = 'C' + RC; }
			else if(SEQ.charAt(x) == 'C') { RC = 'G' + RC; }
			else { RC = 'N' + RC; }
		}
		return RC;
	}
	
	//returns line containing string terminator characters intact
	public static String readLine(BufferedReader READER) throws IOException {
		StringBuilder string = new StringBuilder();
		boolean removeENDLINE = false;
		int x;
		
		while((x = READER.read()) >= 0) {
			if (x == '\n') {
				if(removeENDLINE) {
					removeENDLINE = false;  
				} else {
					string.append((char)'\n');
					break;
				}
			} else {
				removeENDLINE = false;
				string.append((char)x);  
				if (x == '\r') {
					removeENDLINE = true;
					break;
				}
			}
		}
		return string.toString();
	}

	public static void loadConfig(String[] command){
		for (int i = 0; i < command.length; i++) {
			switch (command[i].charAt((1))) {
				case 'b':
					bamFile = new File(command[i + 1]);
					i++;
					break;
				case 'g':
					genome = new File(command[i + 1]);
					i++;
					break;
				case 'o':
					output = new File(command[i + 1]);
					i++;
					break;
				case 'n':
					SEQ = command[i + 1].toUpperCase();
					i++;
					break;
			}
		}
		if(bamFile == null) {
			System.err.println("Invalid BAM File!!!\n");
			printUsage();
			System.exit(1);
		}
		if(genome == null) {
			System.err.println("Invalid Genome File!!!\n");
			printUsage();
			System.exit(1);
		}
		if(!SEQ.equals("A") && !SEQ.equals("T") && !SEQ.equals("G") && !SEQ.equals("C")) {
			System.err.println("Invalid Nucleotide!!!\n");
			printUsage();
			System.exit(1);
		}	
		
		if(output == null) {
			String[] NAME = bamFile.getName().split("\\.");
			output = new File(NAME[0] + "_PSfilter.bam");
		}
				
		System.out.println("-----------------------------------------\nCommand Line Arguments:");
		System.out.println("BAM file: " + bamFile);
		System.out.println("Genome file: " + genome);
		System.out.println("Nucleotide to filter by: " + SEQ);
		System.out.println("Output: " + output);
	}
	
	public static void printUsage() {
		System.err.println("Usage: java -jar filterforPIPseq.jar -b [BAMFile] -g [Genome_Fasta_File] -n [Nucleotide] -o [OutputFile]");
		System.err.println("-----------------------------------------");
		System.err.println("Required: BAM file must be sorted and BAI index must be in same folder as BAM file.");
		System.err.println("\nRequired Parameter:");
		System.err.println("BAM File:\t\t-b\t\tBAM file");
		System.err.println("Genome File:\t\t-g\t\tGenome FASTA File");
		System.err.println("Nucleotide:\t\t-n\t\tNucleotide to filter 5' to the 5' end of the read (A, T, G, or C)");
		
		System.err.println("\nOptional Parameters:");
		System.err.println("Output File Name:\t-o\t\tOutput file");
	}
	
	private static String getTimeStamp() {
		Date date= new Date();
		String time = new Timestamp(date.getTime()).toString();
		return time;
	}
}
