/*
# CLIP-seq CLuster Detection - Tool to detect clusters of an experimental data set of RBP obtained by CLIP-seq protocol.
#
# Created by Paula H. Reyes-Herrera PhD. and Msc. Carlos Andres Sierra on August 2014.
# Copyright (c) 2014 Paula H. Reyes-Herrera PhD. and Msc. Carlos Andres Sierra. Universidad Antonio Narino. All rights reserved.
#
# This file is part of CLIP-seq Cluster Detection.
#
# CLIP-seq Cluster Detection is free software: you can redistribute it and/or modify it under the terms of the 
# GNU General Public License as published by the Free Software Foundation, version 2.
*/


package filesProcess;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Random;
import java.util.Vector;
import java.util.logging.Level;
import java.util.logging.Logger;

import miscellaneous.Functions;
import structures.Background;
import structures.Chromosome;
import structures.Cluster;


/**
 * This class take a sequences-reads set (in SAM Format) to detect and extract the clusters 
 * @author Eng. Paula Reyes, Ph.D. -- M.Sc. Eng. Carlos Sierra
 * Universidad Antonio Narino 
 */
public class CD_FromSam extends Thread 
{

    //Attributes for the class
    private Vector<Chromosome> chromosomes; //Chromosomes information
    private Vector<Cluster> ranking = new Vector<>(); //Ranking of candidate clusters to search motif
    private Vector<Background> nonReads = new Vector<Background>(); //Vector of background
    private String experiment = ""; //DataSet used, extracted from dataset filename
    private String folderResults; //Path for save results
    private boolean bed; //If save clusters in BED format
    private boolean sam; //If save clusters in SAM format
    private int bg_score;
    private int counterSNP = 0; //Amount of SNPs that binding with clusters
    private int counterBackground = 0;
    
    //Ranking parameters
    private int minScore = 0;
    private double minDensity = 0.0;
    
    
    /**
     * Zero-parameters constructor
     */
    public CD_FromSam()
    {
        //Initialize dynamic vector
        this.chromosomes = new Vector<Chromosome>(); //Chromosomes vector
    }
    
    
    /**
     * This method read the file given, extract its information and process data.
     * @param dataset file
     * @param folder of results
     * @param minimum Length-Sequence filter
     * @param minimum Sequences filter
     * @param fileLog
     * @param FASTAs Path
     * @param SNPs Path
     * @param if entropy
     * @param if export to sam
     * @param if export to bed
     * @param Background Path
     * @param if mutations
     * @param if deletions
     * @param if insertions
     * @param window (size of motif)
     * @param amount (ranking clusters)
     * @param setup (setup file)
     * @throws FileNotFoundException
     * @throws InterruptedException
     */
    public void readFile(String file, String folderResults, int minLength, int minSequences, String fileLog, String fastaPath, String snpPath, boolean entropy, boolean sam, boolean bed, String noisePath, String[] mutations, String[] deletions, String[] insertions, int window, int amount, String setup, int bg_score) throws FileNotFoundException, InterruptedException
    {
        this.experiment = file.substring(file.lastIndexOf("/") + 1, file.length()); //Name of data set
        this.folderResults = folderResults; //Folder to save results
        this.bed = bed; //If BED
        this.sam = sam; //If SAM
        this.bg_score = bg_score;
        long startTime = System.nanoTime(); //Start time counter 
        

        //Method Variables  
        String line = "";  //Save the information of one line of file 
        String log = "";  //Save the information of execution 
        
        
        //Read-sequence information
        int start = 0; //Start position of read
        int end = 0; //End position of read
        int seqLenght = 0; //Length of read
        String cigar = ""; //CIGAR of read
        String strand = "*"; //Strand of read
        String sequence = ""; //Sequence of read
        String quality = ""; //Quality of read
        int qmap = 0; //Quality level of mapped process
        int numberMutations = 0; //Number of mutations in read
        String chromosome = "Undef"; //Chromosome of read
        
        
        //Experimental data set information (Use for statistics analysis)
        int contNoValids = 0; //Non-valid reads
        int contValids = 0; //Valid reads
        int contTotals = 0; //Total of reads in file
        int contNoEntrophies = 0; //Number of reads with wrong entropy
        int sequences_fl = 0; //Number of sequences with length minus to filter given. Don't accomplished the filter
        int sequences_fl_vr = 0; //Number of sequences in valid-reads with length minus to filter given
        int sequences_bg = 0; //Number of sequences find in background
        int[] occurrences = new int[4]; //Occurrences of the strands

        
        System.out.println("=====Welcome to CLIP-seq Cluster Detection===== \n\nProcessing...");   log += "Processing...\n";
        
        //Load background
        if(!noisePath.toLowerCase().equals("none"))
        {
        	this.nonReads(noisePath); 
        }
        
        
        //-- Process file --//
        try 
        {
            //Load the file in a buffer
            BufferedReader br = new BufferedReader(new FileReader(file));
            
            //LOG: Describe initial steps
            System.out.println("Separate lines...");        log += "Separate lines...\n";
            System.out.println("Start..."); 	            log += "Start...\n";
            
            //Read first line
            try 
            {
                line = br.readLine();
            } 
            catch (IOException e) 
            {
                System.out.println("ERROR: ClusterDetection.ReadFile(). Can't read first line.");
                log += "ERROR: ClusterDetection.ReadFile(). Can't read first line.\n";
            }
                         
            //LOG
            System.out.println("Header...");      		log += "Header...\n";
                         
            //Read header
            while(line != null && line.substring(0, 3).equals("@HD"))
            {
                try 
                {
                    line = br.readLine();
                } 
                catch (IOException e) 
                {
                    System.out.println("ERROR: ClusterDetection.ReadFile(). Can't read header line.");
                    log += "ERROR: ClusterDetection.ReadFile(). Can't read header line.\n";
                }
            }
             
            
            //LOG: Get chromosomes information
            System.out.println("Chromosomes Information...");		log += "Chromosomes Information...\n";
            
            //Read chromosomes info
            String[] chrInfo; //Array of chomosome info (name, length)
            while(line != null && line.substring(0, 3).equals("@SQ"))
            {
                //Chromosome   Length
                chrInfo = line.split("\t");
               
                if(!chrInfo[1].split(":")[1].substring(0, 4).equals("chrU")) //Avoid chrU
                	this.chromosomes.add(new Chromosome(chrInfo[1].split(":")[1], Integer.parseInt(chrInfo[2].split(":")[1]), folderResults, minSequences, minLength, fastaPath, snpPath));
                
                try 
                {
                    line = br.readLine();
                } 
                catch (IOException e) 
                {
                    System.out.println("ERROR: ClusterDetection.ReadFile(). Can't read chromosomes information.");
                    log += "ERROR: ClusterDetection.ReadFile(). Can't read chromosomes information.\n";
                }
            }
             
            
            //LOG: Process file information
            System.out.println("File Info processing...");
            log += "File Info processing...\n";
                         
            //Read last line before the "reads"
            while(line != null && line.substring(0, 3).equals("@PG"))
            {
                try 
                {
                    line = br.readLine();
                    
                } 
                catch (IOException e) 
                {
                    System.out.println("ERROR: ClusterDetection.ReadFile(). Can't read additional info.");
                    log += "ERROR: ClusterDetection.ReadFile(). Can't read additional info.\n";
                }
            }
            
            
            //When there're no info of chromosomes on SAM format, load some default data
            if(this.chromosomes.size() == 0)
            {
            	this.createChromosomes(folderResults, minSequences, minLength, fastaPath, snpPath);
            }
            
            System.gc(); // Call Garbage Collector
            
            
            
            //Next step: Process reads in file
            System.out.println("\nProcessing reads...");
            log += "\nProcesing reads...\n";
                        
            boolean processRead; //If read can be successful processed
            
            while (line != null) //Line is not empty. We assume that last line in file is empty.
            {
                // ExperimentID   Strand   Chromosome  Start  QMap  CIGAR  -  -  -  Sequence Quality  -  -  NumberMutations
                String[] currentLine = line.split("\t"); //Split the read (line in columns) separated by tabs
                contTotals++; //Increase the amount of reads in file
                
                //Define strand type
                if(currentLine[1].equals("0") || currentLine[1].equals("256"))
                {
                    occurrences[0]++;
                    strand = "+";
                }
                else
                    if(currentLine[1].equals("4"))
                    {
                        occurrences[1]++;
                        strand = "*";
                    }
                    else
                        if(currentLine[1].equals("16") || currentLine[1].equals("272"))
                        {
                            occurrences[2]++;
                            strand = "-";
                        }
                        else
                        {
                            occurrences[3]++;
                            strand = "_";
                        }    
                
                seqLenght = currentLine[9].length(); //Length of the sequence
                
                //Sequences smaller that minLength
                if(seqLenght < minLength)
                {
                	sequences_fl++;
                }
                
                chromosome = currentLine[2]; //Get chromosome of current read
                               
                
                if(!currentLine[2].equals("*") && !currentLine[2].equals("_") && !chromosome.substring(0, 4).equals("chrU")) //Valid if the read is aligned to a chromosome
                {
                	sequence = currentLine[9]; //Sequence of current read
                	processRead = true; //Read can be processed
                    
                    if(entropy) //Define by parameters: Use entropy filter
                    {
                    	//TODO explain constants
                    	processRead = (this.mononucleotideEntrophy(sequence) >= 0.7 && this.dinucleotideEntrophy(sequence) >= 1.5) ? true : false;
                    	
                    	if(!processRead)
                    	{
                    		contNoEntrophies++; //Increase reads discarded by entropy
                    	}
                    }
                    
                    start = Integer.parseInt(currentLine[3]); //Parse value of start position of current read
                    end = start + (seqLenght - 1); //Define end-position of read
                    
                    
                    if(this.nonReads.size() > 0 && processRead) //If there're reads in background, compare an discard
                    {
                    	processRead = this.validateNoiseRead(chromosome, strand, start, end);
                    }	
                    
                    
                    if(processRead) //If read passes filters
                    {
                        contValids++; //Increase valid reads

                        //Valid-sequences smaller that minLength, i.e., not accomplished both length and sequences filters.
                        if(seqLenght < minLength)
                        {
                        	sequences_fl_vr++;
                        }
                        
                        quality = currentLine[10]; //Quality of current read
                        qmap = Integer.parseInt(currentLine[4]); //Quality of mapped of current read
                        
                        //Number of mutations of current read, flexible for TopHap, Bowtie, or anyone
                        for(int i = 11; i < currentLine.length; i++)
                        {
                        	if(currentLine[i].split(":")[0].equals("NM")) //NM is the indicator for this parameter
                        	{
                        		numberMutations = Integer.parseInt(currentLine[i].split(":")[2]);
                        		break;
                        	}
                        }
                         
                        
                        cigar = currentLine[5]; //CIGAR of current read

                        //Move through defined chromosomes
                        for(int i = 0; i < this.chromosomes.size(); i++)
                        {
                            //Find the respective chromosome
                            if(this.chromosomes.get(i).getChromosomeName().equals(chromosome))
                            {
                                //If read has a right strand
                                if(strand.equals("+") || strand.equals("-"))  
                                {
                                	int length = 0; //Length additional by deleted bases
                                	int index = 0; //Index inside sequence
                                	
                                	for(int k = 0; k < cigar.length(); k++)
                    	            {
                    	                if(index == cigar.length()) //If process al CIGAR info, break
                    	                {
                    	                	break;
                    	                }
                    	                
                    	                if(cigar.charAt(k) == 'M') //MATCH not increase length
                    	                {
                    	                	index = k + 1; //Move to index
                    	                }
                    	                else
                    	                    if(cigar.charAt(k) == 'I') //INSERT not increase length
                    	                    {
                    	                        index = k + 1; //Move to index
                    	                    }
                    	                    else
                    	                        if(cigar.charAt(k) == 'D') //DELETION increase length
                    	                        {
                    	                            length += Integer.parseInt(cigar.substring(index, k)); //Add deleted bases
                    	                            index = k + 1; //Move to index
                    	                        }
                    	            }
                                	
                                	//Add read-sequence to respective chromosome
                                	this.chromosomes.get(i).addSequence(currentLine[0], chromosome, start, end + length, strand, qmap, numberMutations, sequence, quality, cigar);
                                }
                            }
                        }
                    }
                    else
                    {
                        sequences_bg++; //Increase reads discarded by background 
                    }
                }
                else
                {
                    contNoValids++; //Increase amount of no-valid reads
                }
                
                
                //Check if still there're non-process lines 
                try 
                {
                    line = br.readLine();
                    if(line == null)
                    {
                    	this.counterBackground++;
                    	
                    	//Print and save values
                        System.out.println("Total reads: " + contTotals + "\tNo valid reads: " + contNoValids + "\tValid reads: " + contValids);
                        System.out.println("Reads Lengths-Filter: " + sequences_fl + "\tValid-reads length filter: " + sequences_fl_vr);
                        System.out.println("Reads Background-Filter: " + this.counterBackground);
                        log += "\nTotal reads: " + contTotals + "\tNo valid reads: " + contNoValids + "\tValid reads: " + contValids; 
                        log += "\nReads Lengths-Filter: " + sequences_fl + "\tValid-reads length filter: " + sequences_fl_vr;
                        log += "\nReads Background-Filter: " + this.counterBackground + "\n";
                        
                        if(entropy) //Add entropies statistics
                        {
                        	System.out.println("Entrophy filter: " + contNoEntrophies);
                        	log += "\tEntrophy filter: " + contNoEntrophies; 
                        }
                        		 
                        log += "\n";
                        
                        br.close();
                        break;
                    }
                } 
                catch (IOException e) 
                {
                	System.out.println("ERROR: ClusterDetection.ReadFile(). Can't save statistical info.");
                    log += "ERROR: ClusterDetection.ReadFile(). Can't save statistical info.\n";
                }
            }
        } 
        finally 
        {
            System.out.println("Closing dataset's file without errors...");
            log += "Closing dataset's without errors...\n";
        }
	
        this.nonReads = null;
        System.gc(); //Call to Garbage Collector
        
        //Get execution time until this point
        double elapsedTimeInSec;
        elapsedTimeInSec = (System.nanoTime() - startTime) * 1.0e-9;
        System.out.println("\n\nElapsed Time: " + elapsedTimeInSec + " seconds\n\nStart Clusters Detection");
        log += "\n\nElapsed Time: " + elapsedTimeInSec + " seconds\n\nProfiling clusters...";
        
        
        //Process choromosomes in parallel
        for(int l = 0; l < this.chromosomes.size(); l++)
        {
            this.chromosomes.get(l).start(); //Start clusters detection for each chromosome
            Thread.sleep(500); //Delay to start next chromosome clusters-detection. Try to avoid concurrence problems.
        }    
        
        int tempScore;
        double tempDensity;
        
        for(int l = 0; l < this.chromosomes.size(); l++)
        {
            this.chromosomes.get(l).join(); //Wait to all chromosomes finish clusters detection
            
            tempScore = this.chromosomes.get(l).getClustersFilter().getMaxScore();
            tempDensity = this.chromosomes.get(l).getClustersFilter().getMaxDensity();
            
            this.minScore = tempScore > this.minScore ? tempScore : this.minScore; 
            this.minDensity = tempDensity > this.minDensity ? tempDensity : this.minDensity;
        } 
        
        this.minScore *= 0.75;
        this.minDensity *= 0.75;
        
        System.out.println("Score: " + this.minScore + "\tDesnsity:" + this.minDensity + "\n");
        
        //Get execution time until this point
        elapsedTimeInSec = (System.nanoTime() - startTime) * 1.0e-9;
        System.out.println("\nElapsed time in process: " + elapsedTimeInSec + " seconds");
    	System.out.println("Saving clusters");
    	log += "\n\nElapsed time in process: " + elapsedTimeInSec + " seconds\n\nSaving clusters";
    	
    	//Start thread that saves general file clusters
       	start();
    	join();
    	
    	
        //Print statistics of execution
        System.out.println("\n\nOccurrences Strand: \n + \t" + occurrences[0] + "\n ind \t" + occurrences[1] + "\n - \t" + occurrences[2] + "\n other \t" + occurrences[3]);
        log += "\n\n\nOccurrences Strand: \n + \t" + occurrences[0] + "\n ind \t" + occurrences[1] + "\n - \t" + occurrences[2] + "\n other \t" + occurrences[3];
        
        System.out.println("\nSNP Clusters: " + this.counterSNP + " clusters");
        log += "\nSNP Clusters: " + this.counterSNP + " clusters\n";
        
        System.out.println("\nBackground: " + sequences_bg + " reads");
        log += "\nSNP Clusters: " + sequences_bg + " reads\n";
        
        //Get execution time until this point
        elapsedTimeInSec = (System.nanoTime() - startTime) * 1.0e-9;
        System.out.println("\n\nElapsed Total Time: " + elapsedTimeInSec + " seconds");
        log += "\n\nElapsed Total Time: " + elapsedTimeInSec + " seconds";
        
        
        
        
        
        //Get trials to make cross-validation: 10 subsets for training, 10 subsets for testing
        // 
        this.getTrials(window);
        double tempTime = elapsedTimeInSec;
        elapsedTimeInSec = (System.nanoTime() - startTime) * 1.0e-9;
        System.out.println("\n\nTrials Building Time: " + (elapsedTimeInSec - tempTime) + " seconds");
        log += "\n\nTrials Building Time: " + (elapsedTimeInSec - tempTime) + " seconds";
        
        //Save log data in a file
        Functions.saveInFile(fileLog, log);
        
        //Stops this thread - Ends of cluster detection process
        this.interrupt();
        
        //Start priority to motif search
        //Priority priority = new Priority(setup);
    }
    
    
    /**
     * 
     * @param folderResults
     * @param minSequences
     * @param minLength
     * @param fastaPath
     * @param snpPath
     */
    public void createChromosomes(String folderResults, int minSequences, int minLength, String fastaPath, String snpPath)
    {
    	this.chromosomes.add(new Chromosome("chr1", 249250621, folderResults, minSequences, minLength, fastaPath, snpPath));
    	this.chromosomes.add(new Chromosome("chr2", 243199373, folderResults, minSequences, minLength, fastaPath, snpPath));
    	this.chromosomes.add(new Chromosome("chr3", 198022430, folderResults, minSequences, minLength, fastaPath, snpPath));
    	this.chromosomes.add(new Chromosome("chr4", 191154276, folderResults, minSequences, minLength, fastaPath, snpPath));
    	this.chromosomes.add(new Chromosome("chr5", 180915260, folderResults, minSequences, minLength, fastaPath, snpPath));
    	this.chromosomes.add(new Chromosome("chr6", 171115067, folderResults, minSequences, minLength, fastaPath, snpPath));
    	this.chromosomes.add(new Chromosome("chr7", 159138663, folderResults, minSequences, minLength, fastaPath, snpPath));
    	this.chromosomes.add(new Chromosome("chr8", 146364022, folderResults, minSequences, minLength, fastaPath, snpPath));
    	this.chromosomes.add(new Chromosome("chr9", 141213431, folderResults, minSequences, minLength, fastaPath, snpPath));
    	this.chromosomes.add(new Chromosome("chr10", 135534747, folderResults, minSequences, minLength, fastaPath, snpPath));
    	this.chromosomes.add(new Chromosome("chr11", 135006516, folderResults, minSequences, minLength, fastaPath, snpPath));
    	this.chromosomes.add(new Chromosome("chr12", 133851895, folderResults, minSequences, minLength, fastaPath, snpPath));
    	this.chromosomes.add(new Chromosome("chr13", 115169878, folderResults, minSequences, minLength, fastaPath, snpPath));
    	this.chromosomes.add(new Chromosome("chr14", 107349540, folderResults, minSequences, minLength, fastaPath, snpPath));
    	this.chromosomes.add(new Chromosome("chr15", 102531392, folderResults, minSequences, minLength, fastaPath, snpPath));
    	this.chromosomes.add(new Chromosome("chr16", 90354753, folderResults, minSequences, minLength, fastaPath, snpPath));
    	this.chromosomes.add(new Chromosome("chr17", 81195210, folderResults, minSequences, minLength, fastaPath, snpPath));
    	this.chromosomes.add(new Chromosome("chr18", 78077248, folderResults, minSequences, minLength, fastaPath, snpPath));
    	this.chromosomes.add(new Chromosome("chr19", 59128983, folderResults, minSequences, minLength, fastaPath, snpPath));
    	this.chromosomes.add(new Chromosome("chr20", 63025520, folderResults, minSequences, minLength, fastaPath, snpPath));
    	this.chromosomes.add(new Chromosome("chr21", 48129895, folderResults, minSequences, minLength, fastaPath, snpPath));
    	this.chromosomes.add(new Chromosome("chr22", 51304566, folderResults, minSequences, minLength, fastaPath, snpPath));
    	this.chromosomes.add(new Chromosome("chrX", 155270560, folderResults, minSequences, minLength, fastaPath, snpPath));
    	this.chromosomes.add(new Chromosome("chrY", 59373566, folderResults, minSequences, minLength, fastaPath, snpPath));
    	this.chromosomes.add(new Chromosome("chrZ", 16571, folderResults, minSequences, minLength, fastaPath, snpPath));
    }
    
    
    /**
     * 
     * @param noisePath
     */
    public void nonReads(String noisePath)
    {
    	BufferedReader br;
        String line;
        String[] currentLine; 
        
        int start;
        int end;
        int score;
        String chromosome;
        String strand;
        
        try 
        {
    		//Load the file in a buffer
    		br = new BufferedReader(new FileReader(noisePath));
            line = br.readLine();
            
            //Read lines
            while(line != null)
            {
            	// ExperimentID   Strand   Chromosome  Start  QMap  CIGAR  -  -  -  Sequence Quality  -  -  NumberMutations
                currentLine = line.split("\t"); //Split the read (line in columns) separate by tabs
            	
                chromosome = currentLine[0];  //Parse the value of start position of current read
                strand = currentLine[5]; //Parse the value of start position of current read
                start = Integer.parseInt(currentLine[1]); //Parse the value of start position of current read
                end = Integer.parseInt(currentLine[2]); //Parse the value of start position of current read
                score = Integer.parseInt(currentLine[4]); //Parse the value of start position of current read
                
                int i = 0;
                
                if(score >= this.bg_score)
                {
	            	for(i = 0; i < this.nonReads.size(); i++)
	                {
	                	if(this.nonReads.get(i).getName().equals(chromosome))
	                	{
	                		this.nonReads.get(i).insert(start, end, strand, score);
	                		break;
	                	}
	                }
	                
	                if(i == this.nonReads.size())
	                {
	                	this.nonReads.add(new Background(chromosome));
	                	this.nonReads.get(i).insert(start, end, strand, score);
	                }
                }
                
                line = br.readLine();
            }
        } 
        catch (IOException e) 
        {
            System.out.println("ERROR: ClusterDetection.ReadNoiseFile(). Can't read first line.");
        }
    }
    
    
    /**
     * 
     * @param chromosome
     * @param strand
     * @param start
     * @param end
     * @return
     */
    public boolean validateNoiseRead(String chromosome, String strand, int start, int end)
    {
    	boolean success = true;
    	
    	for(int i = 0; i < this.nonReads.size(); i++)
    	{
	    	if(this.nonReads.get(i).getName().equals(chromosome))
	    	{
	    		success = !this.nonReads.get(i).search(start, end, strand);
	    		
	    		if(!success)
	    		{
	    			this.counterBackground++;
	    		}
	    		
	    		break;
	    	}
    	}	
    	
    	return success;
    }
    
    
    /**
     * 
     * @param sequence
     * @return 
     */
    public double mononucleotideEntrophy(String sequence)
    {
        double[] probability = new double[4];
        double entrophy = 0.0;
        
        probability[0] = (float)numberOccurrencesMononucleotide(sequence.toUpperCase(), 'A') / (float)sequence.length(); 
        probability[1] = (float)numberOccurrencesMononucleotide(sequence.toUpperCase(), 'C') / (float)sequence.length();    
        probability[2] = (float)numberOccurrencesMononucleotide(sequence.toUpperCase(), 'G') / (float)sequence.length();    
        probability[3] = (float)numberOccurrencesMononucleotide(sequence.toUpperCase(), 'T') / (float)sequence.length();    
        
        
        for(int i = 0; i < 4; i++)
            if(probability[i] != 0.0)
                entrophy += (double)(probability[i] * (double)(Math.log(probability[i])/Math.log(2)+1e-10));
        
        
        entrophy *= (-1.0);
      
        return entrophy;
    }
 
    
    /**
     * 
     * @param sequence
     * @param base
     * @return 
     */
    public int numberOccurrencesMononucleotide(String sequence, char base)
    {
        int counter = 0;
        
        for(int i = 0; i < sequence.length(); i++)
           if(sequence.charAt(i) == base)
               counter++;
        
        return counter;
    }
    
    
    /**
     * 
     * @param sequence
     * @return 
     */
    public double dinucleotideEntrophy(String sequence)
    {
        double entrophy = 0.0;
        double[] probability = new double[12];
        
        probability[0] = (double)numberOccurrencesDinucleotide(sequence.toUpperCase(), "AC") / (double)sequence.length();    
        probability[1] = (double)numberOccurrencesDinucleotide(sequence.toUpperCase(), "AG") / (double)sequence.length();    
        probability[2] = (double)numberOccurrencesDinucleotide(sequence.toUpperCase(), "AT") / (double)sequence.length();    
        probability[3] = (double)numberOccurrencesDinucleotide(sequence.toUpperCase(), "CA") / (double)sequence.length();    
        probability[4] = (double)numberOccurrencesDinucleotide(sequence.toUpperCase(), "CG") / (double)sequence.length();    
        probability[5] = (double)numberOccurrencesDinucleotide(sequence.toUpperCase(), "CT") / (double)sequence.length();    
        probability[6] = (double)numberOccurrencesDinucleotide(sequence.toUpperCase(), "GA") / (double)sequence.length();    
        probability[7] = (double)numberOccurrencesDinucleotide(sequence.toUpperCase(), "GC") / (double)sequence.length();    
        probability[8] = (double)numberOccurrencesDinucleotide(sequence.toUpperCase(), "GT") / (double)sequence.length();    
        probability[9] = (double)numberOccurrencesDinucleotide(sequence.toUpperCase(), "TA") / (double)sequence.length();    
        probability[10] = (double)numberOccurrencesDinucleotide(sequence.toUpperCase(), "TC") / (double)sequence.length();    
        probability[11] = (double)numberOccurrencesDinucleotide(sequence.toUpperCase(), "TG") / (double)sequence.length();    
        
        for(int i = 0; i < 12; i++)
            if(probability[i] != 0.0)
                entrophy += (double)(probability[i] * (double)(Math.log(probability[i])/Math.log(2)+1e-10));
        
        entrophy *= (-1);
        
        return entrophy;
    }
    
    
    /**
     * 
     * @param sequence
     * @param base
     * @return 
     */
    public int numberOccurrencesDinucleotide(String sequence, String base)
    {
        int counter = 0;
        
        for(int i = 0; i < (sequence.length() - 1); i++)
           if(sequence.substring(i, i+2).compareTo(base) == 0)
               counter++;
        
        return counter;
    }
 
    
    
    /**
     * 
     */
    public void randomDouble()
    {
    	Cluster swap = null;
    	Random rd = new Random();
    	int tempIndex;
    	int n = this.ranking.size();
    	
    	for(int i = 0; i < n; i++)
    	{
    		tempIndex = rd.nextInt(n);
    		
    		//Make swap
    		swap = this.ranking.get(i);
    		this.ranking.get(i).change(this.ranking.get(tempIndex)); 
    		this.ranking.get(tempIndex).change(swap);
    	}
    	
    	for(int i = (n - 1); i >= 0; i--)
    	{
    		tempIndex = rd.nextInt(n);
    		
    		//Make swap
    		swap = this.ranking.get(i);
    		this.ranking.get(i).change(this.ranking.get(tempIndex)); 
    		this.ranking.get(tempIndex).change(swap);
    	}
    }
    
    
    
    /**
     * This method is used to save the clusters information in a file with UAN format.
     */
    public void saveTxt()
    {
    	int counterReads = 1;
    	int counterReadsFilter = 1;
    	String fileTotalTxt = folderResults + "/Clusters_total.txt";
    	String fileTotalSNP = folderResults + "/Clusters_total_SNPs.txt";
    	String fileFilterTxt = folderResults + "/Clusters_filtered.txt";
    	String fileFilterSNP = folderResults + "/Clusters_filtered_SNPs.txt";
    	
    	
    	//-------------------- Total Files --------------------//
    	File fileSaveTotalTxt = new File(fileTotalTxt);
    	FileWriter fwTotalTxt;
        if (!fileSaveTotalTxt.exists()) 
        {
            try 
            {
                fileSaveTotalTxt.createNewFile();
            } 
            catch (IOException ex) 
            {
                Logger.getLogger(Chromosome.class.getName()).log(Level.SEVERE, null, ex);
                System.out.println("\n\nFile " + fileTotalTxt);
                System.exit(0);
            }
        }
        
        File fileSaveTotalSNP = new File(fileTotalSNP);
    	FileWriter fwTotalSNP;
        if (!fileSaveTotalSNP.exists()) 
        {
            try 
            {
                fileSaveTotalSNP.createNewFile();
            } 
            catch (IOException ex) 
            {
                Logger.getLogger(Chromosome.class.getName()).log(Level.SEVERE, null, ex);
                System.out.println("\n\nFile " + fileTotalSNP);
                System.exit(0);
            }
        }
        
        
      //-------------------- Filter Files --------------------//
    	File fileSaveFilterTxt = new File(fileFilterTxt);
    	FileWriter fwFilterTxt;
        if (!fileSaveFilterTxt.exists()) 
        {
            try 
            {
                fileSaveFilterTxt.createNewFile();
            } 
            catch (IOException ex) 
            {
                Logger.getLogger(Chromosome.class.getName()).log(Level.SEVERE, null, ex);
                System.out.println("\n\nFile " + fileFilterTxt);
                System.exit(0);
            }
        }
        
        File fileSaveFilterSNP = new File(fileFilterSNP);
    	FileWriter fwFilterSNP;
        if (!fileSaveFilterSNP.exists()) 
        {
            try 
            {
                fileSaveFilterSNP.createNewFile();
            } 
            catch (IOException ex) 
            {
                Logger.getLogger(Chromosome.class.getName()).log(Level.SEVERE, null, ex);
                System.out.println("\n\nFile " + fileFilterSNP);
                System.exit(0);
            }
        }
        
        
        
        try 
        {
        	fwTotalTxt = new FileWriter(fileSaveTotalTxt.getAbsoluteFile(), true);
			BufferedWriter bwTotalTxt = new BufferedWriter(fwTotalTxt);
        	
        	fwFilterTxt = new FileWriter(fileSaveFilterTxt.getAbsoluteFile(), true);
			BufferedWriter bwFilterTxt = new BufferedWriter(fwFilterTxt);
			
			fwTotalSNP = new FileWriter(fileSaveTotalSNP.getAbsoluteFile(), true);
			BufferedWriter bwTotalSNP = new BufferedWriter(fwTotalSNP);
			
			fwFilterSNP = new FileWriter(fileSaveFilterSNP.getAbsoluteFile(), true);
			BufferedWriter bwFilterSNP = new BufferedWriter(fwFilterSNP);
		    
	    	for(int l = 0; l < this.chromosomes.size(); l++)
	        {    
	    		for(int n = 0; n < this.chromosomes.get(l).getClusters().getSet().size(); n++)
	            {
	    			try 
	    			{
						bwTotalTxt.write(counterReads + "\t" + (n + 1) + "\t" + this.chromosomes.get(l).getClusters().getSet().get(n).toString(this.experiment + "." + counterReads));
						
						if(this.chromosomes.get(l).getClusters().getSet().get(n).SNPcontact())
						{
							bwTotalSNP.write(counterReads + "\t" + (n + 1) + "\t" + this.chromosomes.get(l).getClusters().getSet().get(n).toString(this.experiment + "." + counterReads));
						}
						else
						{
							counterSNP++;
						}
					} 
	    			catch (IOException e) 
	    			{
						e.printStackTrace();
					}
	            
	                counterReads++;
	            }
	    		
	    		for(int n = 0; n < this.chromosomes.get(l).getClustersFilter().getSet().size(); n++)
	            {
	    			try 
	    			{
						bwFilterTxt.write(counterReadsFilter + "\t" + (n + 1) + "\t" + this.chromosomes.get(l).getClustersFilter().getSet().get(n).toString(this.experiment + "." + counterReadsFilter));
						
						if(this.chromosomes.get(l).getClustersFilter().getSet().get(n).SNPcontact())
						{
							bwFilterSNP.write(counterReadsFilter + "\t" + (n + 1) + "\t" + this.chromosomes.get(l).getClustersFilter().getSet().get(n).toString(this.experiment + "." + counterReadsFilter));
						}
					} 
	    			catch (IOException e) 
	    			{
						e.printStackTrace();
					}
	            
	                counterReadsFilter++;
	            }    
	        }
	    	
	    	bwTotalTxt.close();
	    	bwTotalSNP.close();
	    	bwFilterTxt.close();
	    	bwFilterSNP.close();
		} 
        catch (IOException e1) 
        {
			e1.printStackTrace();
		}
    }
    
    
    /**
     * This method is used to save the clusters information in a file with BED format. 
     */
    public void saveBed()
    {
    	int counterReads = 1;
    	int counterReadsFilter = 1;
    	String fileTotal = folderResults + "/Clusters_total.bed";
    	String fileFilter = folderResults + "/Clusters_filtered.bed";
    	
    	File fileSaveTotal = new File(fileTotal);
    	FileWriter fwTotal;
        if (!fileSaveTotal.exists()) 
        {
            try 
            {
                fileSaveTotal.createNewFile();
            } 
            catch (IOException ex) 
            {
                Logger.getLogger(Chromosome.class.getName()).log(Level.SEVERE, null, ex);
                System.out.println("\n\nFile " + fileTotal);
                System.exit(0);
            }
        }
        
        File fileSaveFilter = new File(fileFilter);
    	FileWriter fwFilter;
        if (!fileSaveFilter.exists()) 
        {
            try 
            {
                fileSaveFilter.createNewFile();
            } 
            catch (IOException ex) 
            {
                Logger.getLogger(Chromosome.class.getName()).log(Level.SEVERE, null, ex);
                System.out.println("\n\nFile " + fileFilter);
                System.exit(0);
            }
        }
        
        
        try 
        {
			fwTotal = new FileWriter(fileSaveTotal.getAbsoluteFile(), true);
			BufferedWriter bwTotal = new BufferedWriter(fwTotal);
			
			fwFilter = new FileWriter(fileSaveFilter.getAbsoluteFile(), true);
			BufferedWriter bwFilter = new BufferedWriter(fwFilter);
			
		    
	    	for(int l = 0; l < this.chromosomes.size(); l++)
	        {    
	    		for(int n = 0; n < this.chromosomes.get(l).getClusters().getSet().size(); n++)
	            {
	    			try 
	    			{
						bwTotal.write(this.chromosomes.get(l).getClusters().getSet().get(n).toString_BED(this.experiment + "." + counterReads));
					} 
	    			catch (IOException e) 
	    			{
						e.printStackTrace();
					}
	            
	                counterReads++;
	            }
	    		
	    		for(int n = 0; n < this.chromosomes.get(l).getClustersFilter().getSet().size(); n++)
	            {
	    			try 
	    			{
						bwFilter.write(this.chromosomes.get(l).getClustersFilter().getSet().get(n).toString_BED(this.experiment + "." + counterReadsFilter));
					} 
	    			catch (IOException e) 
	    			{
						e.printStackTrace();
					}
	            
	                counterReadsFilter++;
	            }    
	        }
	    	
	    	bwTotal.close();
	    	bwFilter.close();
		} 
        catch (IOException e1) 
        {
			e1.printStackTrace();
		}
    }
    
    
    /**
     * This method is used to save the clusters information in a file with SAM format.
     */
    public void saveSam()
    {
    	int counterReads = 1;
    	int counterReadsFilter = 1;
    	String fileTotal = folderResults + "/Clusters_total.sam";
    	String fileFilter = folderResults + "/Clusters_filtered.sam";
    	
    	File fileSaveTotal = new File(fileTotal);
    	FileWriter fwTotal;
        if (!fileSaveTotal.exists()) 
        {
            try 
            {
                fileSaveTotal.createNewFile();
            } 
            catch (IOException ex) 
            {
                Logger.getLogger(Chromosome.class.getName()).log(Level.SEVERE, null, ex);
                System.out.println("\n\nFile " + fileTotal);
                System.exit(0);
            }
        }
        
        File fileSaveFilter = new File(fileFilter);
    	FileWriter fwFilter;
        if (!fileSaveFilter.exists()) 
        {
            try 
            {
                fileSaveFilter.createNewFile();
            } 
            catch (IOException ex) 
            {
                Logger.getLogger(Chromosome.class.getName()).log(Level.SEVERE, null, ex);
                System.out.println("\n\nFile " + fileFilter);
                System.exit(0);
            }
        }
        
        
        try 
        {
			fwTotal = new FileWriter(fileSaveTotal.getAbsoluteFile(), true);
			BufferedWriter bwTotal = new BufferedWriter(fwTotal);
			
			fwFilter = new FileWriter(fileSaveFilter.getAbsoluteFile(), true);
			BufferedWriter bwFilter = new BufferedWriter(fwFilter);
		    
	    	for(int l = 0; l < this.chromosomes.size(); l++)
	        {    
	    		for(int n = 0; n < this.chromosomes.get(l).getClusters().getSet().size(); n++)
	            {
	    			try 
	    			{
						bwTotal.write(this.chromosomes.get(l).getClusters().getSet().get(n).toString_SAM(this.experiment + "." + counterReads));
					} 
	    			catch (IOException e) 
	    			{
						e.printStackTrace();
					}
	            
	                counterReads++;
	            }    
	    		
	    		for(int n = 0; n < this.chromosomes.get(l).getClustersFilter().getSet().size(); n++)
	            {
	    			try 
	    			{
						bwFilter.write(this.chromosomes.get(l).getClustersFilter().getSet().get(n).toString_SAM(this.experiment + "." + counterReadsFilter));
					} 
	    			catch (IOException e) 
	    			{
						e.printStackTrace();
					}
	            
	                counterReadsFilter++;
	            }
	        }
	    	
	    	bwTotal.close();
	    	bwFilter.close();
		} 
        catch (IOException e1) 
        {
			e1.printStackTrace();
		}
    }
    
    
    /**
     * This method is used to save the clusters information in a file with SAM format.
     */
    public void saveRanking()
    {
    	String fileRanking = folderResults + "/Ranking.bed";
    	
    	File fileSaveRanking = new File(fileRanking);
    	FileWriter fwRanking;
        if (!fileSaveRanking.exists()) 
        {
            try 
            {
                fileSaveRanking.createNewFile();
            } 
            catch (IOException ex) 
            {
                Logger.getLogger(Chromosome.class.getName()).log(Level.SEVERE, null, ex);
                System.out.println("\n\nFile " + fileRanking);
                System.exit(0);
            }
        }
        
        try 
        {
			fwRanking = new FileWriter(fileSaveRanking.getAbsoluteFile(), true);
			BufferedWriter bwRanking = new BufferedWriter(fwRanking);
			
			
	    	for(int l = 0; l < this.ranking.size(); l++)
	        {    
	    		try 
	    		{
					bwRanking.write(this.ranking.get(l).toString_BED("Ranking." + (l + 1)) + "\n");
				} 
	    		catch (IOException e) 
	    		{
					e.printStackTrace();
				}
	            
	        }
	    	
	    	bwRanking.close();
	    } 
        catch (IOException e1) 
        {
			e1.printStackTrace();
		}
    }
    
    /**
     * 
     * @param window
     */
    public void getTrials(int window)
    {
    	int length = this.ranking.size();
    	int subset = length / 10;
    	
    	String pathFasta = folderResults + "/FASTA/";
    	String pathPrior = folderResults + "/PRIORS/";
    	
    	File fileFasta = new File(pathFasta);
    	fileFasta.mkdir();
    	
    	File filePrior = new File(pathPrior);
    	filePrior.mkdir();
        
    	
    	for(int i = 0; i < 10; i++)
    	{
    		saveTrials(i + 1, i * subset, (i * subset) + (subset - 1), window);
    	}
    	
    	//saveRankings(1, 0, length, window);
    }
    
    
    /**
     * 
     * @param set
     * @param start
     * @param end
     * @param window
     */
    public void saveTrials(int set, int start, int end, int window)
    {
    	String fileFastaTraining = folderResults + "/FASTA/trial_" + set + "_training.fasta";
    	String filePriorTraining = folderResults + "/PRIORS/trial_" + set + "_training.prior";
    	String fileFastaTesting = folderResults + "/FASTA/trial_" + set + "_testing.fasta";
    	String filePriorTesting = folderResults + "/PRIORS/trial_" + set + "_testing.prior";
    	String fileFastaTestingNeg = folderResults + "/FASTA/trial_" + set + "_testing_negative.fasta";
    	
    	File fileSaveFastaTraining = new File(fileFastaTraining);
    	FileWriter fwFastaTraining;
        if (!fileSaveFastaTraining.exists()) 
        {
            try 
            {
                fileSaveFastaTraining.createNewFile();
            } 
            catch (IOException ex) 
            {
                Logger.getLogger(Chromosome.class.getName()).log(Level.SEVERE, null, ex);
                System.out.println("\n\nFile " + fileFastaTraining);
                System.exit(0);
            }
        }
        
        File fileSaveFastaTesting = new File(fileFastaTesting);
    	FileWriter fwFastaTesting;
        if (!fileSaveFastaTesting.exists()) 
        {
            try 
            {
                fileSaveFastaTesting.createNewFile();
            } 
            catch (IOException ex) 
            {
                Logger.getLogger(Chromosome.class.getName()).log(Level.SEVERE, null, ex);
                System.out.println("\n\nFile " + fileFastaTesting);
                System.exit(0);
            }
        }
        
        File fileSaveFastaTestingNeg = new File(fileFastaTestingNeg);
    	FileWriter fwFastaTestingNeg;
        if (!fileSaveFastaTestingNeg.exists()) 
        {
            try 
            {
                fileSaveFastaTestingNeg.createNewFile();
            } 
            catch (IOException ex) 
            {
                Logger.getLogger(Chromosome.class.getName()).log(Level.SEVERE, null, ex);
                System.out.println("\n\nFile " + fileFastaTestingNeg);
                System.exit(0);
            }
        }
        
        
        File fileSavePriorTraining = new File(filePriorTraining);
    	FileWriter fwPriorTraining;
        if (!fileSavePriorTraining.exists()) 
        {
            try 
            {
                fileSavePriorTraining.createNewFile();
            } 
            catch (IOException ex) 
            {
                Logger.getLogger(Chromosome.class.getName()).log(Level.SEVERE, null, ex);
                System.out.println("\n\nFile " + filePriorTraining);
                System.exit(0);
            }
        }
        
        File fileSavePriorTesting = new File(filePriorTesting);
    	FileWriter fwPriorTesting;
        if (!fileSavePriorTesting.exists()) 
        {
            try 
            {
                fileSavePriorTesting.createNewFile();
            } 
            catch (IOException ex) 
            {
                Logger.getLogger(Chromosome.class.getName()).log(Level.SEVERE, null, ex);
                System.out.println("\n\nFile " + filePriorTesting);
                System.exit(0);
            }
        }
        
        
        try 
        {
			fwFastaTraining = new FileWriter(fileSaveFastaTraining.getAbsoluteFile(), true);
			BufferedWriter bwFastaTraining = new BufferedWriter(fwFastaTraining);
			
			fwFastaTesting = new FileWriter(fileSaveFastaTesting.getAbsoluteFile(), true);
			BufferedWriter bwFastaTesting = new BufferedWriter(fwFastaTesting);
			
			fwFastaTestingNeg = new FileWriter(fileSaveFastaTesting.getAbsoluteFile(), true);
			BufferedWriter bwFastaTestingNeg = new BufferedWriter(fwFastaTestingNeg);
			
			fwPriorTraining = new FileWriter(fileSavePriorTraining.getAbsoluteFile(), true);
			BufferedWriter bwPriorTraining = new BufferedWriter(fwPriorTraining);
			
			fwPriorTesting = new FileWriter(fileSavePriorTesting.getAbsoluteFile(), true);
			BufferedWriter bwPriorTesting = new BufferedWriter(fwPriorTesting);
		    
						
			for(int l = 0; l < this.ranking.size(); l++)
	        {    
	    		if(l >= start && l <= end)
	    		{
		    		try 
		    		{
		    			bwFastaTesting.write(this.ranking.get(l).toString_Fasta());
						bwPriorTesting.write(this.ranking.get(l).toString_Prior(window));
						
						//Negative set
						bwFastaTestingNeg.write(getNegative(this.ranking.get(l)));
					} 
		    		catch (IOException e) 
		    		{
						e.printStackTrace();
					}
	    		}
	    		else
	    		{
	    			try 
		    		{
						bwFastaTraining.write(this.ranking.get(l).toString_Fasta());
						bwPriorTraining.write(this.ranking.get(l).toString_Prior(window));
					} 
		    		catch (IOException e) 
		    		{
						e.printStackTrace();
					}
	    		}
	        }    
			
			bwFastaTraining.close();
	    	bwPriorTraining.close();
	    	
	    	bwFastaTesting.close();
	    	bwFastaTestingNeg.close();
	    	bwPriorTesting.close();
		} 
        catch (IOException e1) 
        {
			e1.printStackTrace();
		}
    }

    
    /**
     * 
     * @param cluster
     * @return
     */
    public String getNegative(Cluster cluster)
    {
    	String negative = "";
    	String seq = cluster.getSequence();
    	double[] weights = new double[5];
    	
    	for(int i = 0; i < seq.length(); i++)
    	{
    		switch(seq.charAt(i))
    		{
    			case 'A':	weights[0] += 1.0;	break;
    			case 'C':	weights[1] += 1.0;	break;
    			case 'G':	weights[2] += 1.0;	break;
    			case 'T':	weights[3] += 1.0;	break;
    			case 'N':	weights[4] += 1.0;	break;
    		}
    	}
    	
    	for(int i = 0; i < 5; i++)
    	{
    		weights[i] /= (double) seq.length();
    	}
    	
    	for(int i = 1; i < 5; i++)
    	{
    		weights[i] += weights[i - 1];
    	}
    	
    	double roulette;
    	Random rd = new Random();
    	
    	for(int i = 0; i < seq.length(); i++)
    	{
    		roulette = rd.nextDouble();
    		
    		if(roulette < weights[0])
    		{
    			negative += "A";
    		}
    		else
    		{
    			if(roulette < weights[1])
        		{
        			negative += "C";
        		}
        		else
        		{
        			if(roulette < weights[2])
            		{
            			negative += "G";
            		}
            		else
            		{
            			if(roulette < weights[3])
                		{
                			negative += "T";
                		}
                		else
                		{
                			negative += "N";
                		}
            		}
        		}
    		}
    	}
    	
    	return negative;
    }
    
    
    /**
     * 
     * @param amount
     */
    public void createRanking()
    {
    	int temp = 0;
    	
    	for(int i = 0; i < this.chromosomes.size(); i++)
    	{
    		for(int j = 0; j < this.chromosomes.get(i).getClustersFilter().getSet().size(); j++)
    		{
    			if(this.chromosomes.get(i).getClustersFilter().getSet().get(j).getMutationsRanking() >= this.minScore ||
    				this.chromosomes.get(i).getClustersFilter().getSet().get(j).getDensityRanking() >= this.minDensity)
    			{
    				temp = findIndexRanking(this.chromosomes.get(i).getClustersFilter().getSet().get(j).getMutationsRanking());
    				this.ranking.add(temp, this.chromosomes.get(i).getClustersFilter().getSet().get(j));
    			}
    		}
    	}
    	
    	this.randomDouble();
    	
    	this.saveRanking();
    }
    
    
    public int findIndexRanking(int score)
    {
    	int min = 0; //Current lower limit of search segment
    	int max = this.ranking.size(); //Current upper limit of search segment
    	int size = this.ranking.size(); //Size of the cluster
 		int middle; //Split point of search segment
 		int index = 0;
			
        do //Use binary search to find the respective index
        {
            middle = (min + max) / 2;

            if (middle >= 0 && (middle < size) && (this.ranking.get(middle).getMutationsRanking() < score))
            	min = middle + 1; //Discard the lower part
            else 
            	max = middle; //Discard the upper part
                
            if(min == max)
                index = min;
        }
        while (min < max);
			
        return index;
    }
    
    //Threads
    @Override
    public void run()
    {
    	Thread rank = new Thread ()
    	{
    		@Override
    		public void run()
    		{
    			createRanking();
    		}
    	};
    	
    	rank.start();
    	
    	
    	Thread txt = new Thread ()
    	{
    		@Override
    		public void run()
    		{
    			saveTxt();
    		}
    	};
    	
    	txt.start();
    	                 	
    	/**
    	 * Thread to save in BED format
    	 */
    	Thread bedT = new Thread () 
    	{
	    	@Override
    		public void run()
	    	{
	    		if(bed)
	        	{
	    			saveBed();
	    		}
	    		
	    	}	
	    };
	    	
	    bedT.start();
    	
	    /**
	     * Thread to save in SAM format
	     */
    	Thread samT = new Thread () 
    	{
    		@Override
    		public void run()
    		{
    			if(sam)
    	    	{
    	    		saveSam();
    			}
    		}	
    	};
    		
    	samT.start();
    	
    	
    	/**
    	 * Finalize the threads when these finish their process.
    	 */
    	try 
    	{
			rank.join();
			rank.interrupt();
    		
    		txt.join();
			txt.interrupt();
			
			bedT.join();
			bedT.interrupt();
			
			samT.join();
			samT.interrupt();
		} 
    	catch (InterruptedException e) 
    	{
			e.printStackTrace();
		}
    }
}