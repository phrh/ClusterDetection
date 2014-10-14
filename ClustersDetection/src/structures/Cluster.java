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


package structures;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Vector;



/**
 * This class represents the behavior of a cluster
 * @author Eng. Paula Reyes, Ph.D. -- M.Sc. Eng. Carlos Sierra
 * Antonio Narino University
 */
public class Cluster extends Thread{

    //Attributes
    private Vector<Read> sequences = null; //Sequences inside the cluster
    private int minPosition; // Start position of the cluster
    private int maxPosition; //End position of the cluster
    private String strand; //Strand associated at the cluster
    private String chromosome; //Chromosome associated at the cluster
	
    //General Indicators
    private int minClusterLength = 1000; //Minimum length of a cluster
    private int maxClusterLength = 0; //Maximum length of a cluster
    private int minClusterSequences = 1000; //Minimum number of sequences in a cluster
    private int maxClusterSequences = 0; //Maximum number of sequences in a cluster
    private int numberMutations = 0; //NUmber of Mutations
    private int numberMutationsRevelant = 0; //Number of mutations with SNPs filter
    
    //Profile
    private String sequence = "";  //Cluster sequence
    private int coverage[] = null; //Coverage vector
    private int mutation[] = null; //General mutations vector
    private int undefined[] = null; //General mutations vector
    private int mutations[][] = null; //Mutations vector by mutation (A-C, A-G, A-T, ..., T-G)
    private int insertions[][] = null; //Insertions vector by nucleotide base (A, C, G, T)
    private int deletions[][] = null; //Deletions vector by nucleotide base (A, C, G, T)

    
    private String logSNP = "";
    
    
    /**
     * Zero-parameters constructor
     */
    public Cluster() 
    {
        this.sequences = null;
        this.minPosition = 0;
        this.maxPosition = 0;
        this.chromosome = null;
    }

	
    /**
     * Constructor with parameters
     * @param minPosition
     * @param maxPosition
     * @param strand
     * @param chromosome 
     */
    public Cluster(int minPosition, int maxPosition, String strand, String chromosome)
    {
    	this.minPosition = minPosition;
    	this.maxPosition = maxPosition;
    	this.strand = strand;
        this.chromosome = chromosome;
        this.sequences = new Vector<Read>();
    }

    
    /**
     * @return the sequence
     */
    public String getSequence() 
    {
        return sequence;
    }

    
    /**
     * @param sequence the sequence to set
     */
    public void setSequence(String sequence) 
    {
        this.sequence = sequence;
    }


    /**
     * 
     * @return minimumPosition
     */
    public int getMinPosition() 
    {
    	return minPosition;
    }


    /**
     * @param minPosition the minimumPosition to set
     */
    public void setMinPosition(int minPosition) 
    {
    	this.minPosition = minPosition;
    }


    /**
     * @return the maxPosition
     */
    public int getMaxPosition() 
    {
    	return maxPosition;
    }


    /**
     * @param maxPosition the maxPosition to set
     */
    public void setMaxPosition(int maxPosition) 
    {
    	this.maxPosition = maxPosition;
    }
	
	
    /**
     * @return the strand
     */
    public String getStrand() 
    {
    	return this.strand;
    }


    /**
     * @param strand the strand to set
     */
    public void setStrand(String strand) 
    {
    	this.strand = strand;
    }


       /**
     * @return the chromosome
     */
    public String getChromosome() 
    {
        return chromosome;
    }

    
    /**
     * @param chromosome the chromosome to set
     */
    public void setChromosome(String chromosome) 
    {
        this.chromosome = chromosome;
    }
    
    
    /**
     * @return the minClusterLength
     */
    public int getMinClusterLength() {
        return minClusterLength;
    }

    /**
     * @param minClusterLength the minClusterLength to set
     */
    public void setMinClusterLength(int minClusterLength) 
    {
        this.minClusterLength = minClusterLength;
    }

    
    /**
     * @return the maxClusterLength
     */
    public int getMaxClusterLength() 
    {
        return maxClusterLength;
    }

    
    /**
     * @param maxClusterLength the maxClusterLength to set
     */
    public void setMaxClusterLength(int maxClusterLength) 
    {
        this.maxClusterLength = maxClusterLength;
    }

    
    /**
     * @return the minClusterSequences
     */
    public int getMinClusterSequences() 
    {
        return minClusterSequences;
    }

    
    /**
     * @param minClusterSequences the minClusterSequences to set
     */
    public void setMinClusterSequences(int minClusterSequences) 
    {
        this.minClusterSequences = minClusterSequences;
    }

    
    /**
     * @return the maxClusterSequences
     */
    public int getMaxClusterSequences() 
    {
        return maxClusterSequences;
    }

    
    /**
     * @param maxClusterSequences the maxClusterSequences to set
     */
    public void setMaxClusterSequences(int maxClusterSequences) 
    {
        this.maxClusterSequences = maxClusterSequences;
    }
    
 
    /**
     * 
     * @return 
     */
    public int length()
    {
        //System.out.println("Length: " + (this.maxPosition - this.minPosition + 1));
    	return this.maxPosition - this.minPosition + 1;
    }
    
    
    /**
     * This method take the parameters given (relevant information of a sequence) and determinate 
     * if the sequence is one of the set of sequences of current cluster
     * @param startPosition
     * @param endPosition
     * @param strand
     * @param chromosome
     * @return if the sequence must be in the current cluster 
     */
    public boolean isSequence(int startPosition, int endPosition, String strand, String chromosome)
    {
    	boolean successful = false; //Default is FALSE. We assume that the sequence is not in the cluster.
		
    	
    	
    	if(this.strand.equals(strand) && this.chromosome.equals(chromosome)) //Validate the strand and chromosome
            if(startPosition <= this.minPosition && endPosition >= this.minPosition) //Sequence is at first part of the cluster
                successful = true;
            else
                if(startPosition <= this.maxPosition && endPosition >= this.maxPosition) //Sequence is at last part of the cluster
                    successful = true;
                else
                    if(startPosition >= this.minPosition && endPosition <= this.maxPosition) //Sequence is inside the cluster
                        successful = true;
                    else
                        if(startPosition <= this.maxPosition && endPosition >= this.maxPosition) //Cluster is inside the sequence
                            successful = true;
                        else
                        	if(startPosition == this.minPosition || endPosition == this.maxPosition) //Cluster is inside the sequence
                        	{
                        		successful = true;
                        		
                        		
                        	}
		
    		
    	
        return successful;  
    }
	
    
    /**
     * This method returns all sequences in the cluster.
     * @return sequences
     */
    public Vector<Read> getSequences()
    {
    	try
    	{
            return this.sequences;
    	}
    	catch (Exception e)
    	{
            System.out.println("ERROR: Cluster.getSequences(). Can't return the sequences in the cluster.");
            return null;
    	}
    }
	
	
    /**
     * This method returns one specific sequence in the cluster.
     * @param index
     * @return sequence in the given index
     */
    public Read getSequence(int index)
    {
        try
	{
            return this.sequences.get(index);
        }
	catch (Exception e)
	{
            System.out.println("ERROR: Cluster.getSequences(index). Can't return the specified sequence in the cluster.");
            return null;
        }
    }
	
	
    /**
     * This method add a sequence inside the cluster. 
     * If the sequence yet exists just increase the occurrences amount of that sequence;
     * else, add a new sequence with just one occurrence.
     * @param sequence
     * @return successful notification
     */
    public boolean addSequence(Read sequence)
    {
	boolean succesfully = false;
		
    	try 
    	{
            int index = this.existsSequence(sequence.getStart(), sequence.getEnd(), sequence.getSequence());
			
            if(index == -1)
            {
            	//Add a new sequence
                this.sequences.add(sequence);
				
                //Validate minimum position in the cluster
                this.minPosition = sequence.getStart() < this.minPosition ? sequence.getStart(): this.minPosition;  
				
                //Validate maximum position in the cluster
                this.maxPosition = sequence.getEnd() > this.maxPosition ? sequence.getEnd(): this.maxPosition;
            }
            else
            {
                //Increase occurrence of a existent sequence
                this.sequences.get(index).increaseOccurrence();
            }
			
            succesfully = true;
    	}
    	catch(Exception e)
    	{
            System.out.println("ERROR: Cluster.addSequence(). Can't add sequence at the cluster.");
        }
		
        return succesfully;
    }
	
	
    /**
     * This method verify if the sequence exists in the cluster or not.
     * @param start
     * @param end
     * @return if the sequence is yet in the cluster
     */
    public int existsSequence(int start, int end, String sequence)
    {
    	int index = -1;
		
    	for(int i = 0; i < this.sequences.size(); i++)
    		if(this.sequences.get(i).getStart() == start && this.sequences.get(i).getEnd() == end && this.sequences.get(i).getSequence().equals(sequence))
                index = i;
                
    	return index;
    }
	
    
    /**
     * This method returns the amount of sequences in the cluster, including the occurrences per sequence.
     * @return amount of sequences in cluster
     */
    public int getAmountSequences()
    {
        int cont = 0;
            
        for(int i = 0; i < this.sequences.size(); i++)
            cont += this.sequences.get(i).getOccurrences();
           
        return cont;
    }

    
    /**
     * This method delete a sequence in the given index.
     * @param index
     * @return successful notification
     */
    @SuppressWarnings("unused")
    private boolean deleteSequence(int index)
    {
        try
        {
            this.sequences.remove(index);
            return true;
        }
        catch (Exception e)
        {
            System.out.println("ERROR: Cluster.deleteSequence(). Index isn't the bounds or the sequence can't be removed.");
            return false;
        }
    }

    
    /**
     * @return the coverage
     */
    public int[] getCoverage() 
    {
        return coverage;
    }

    
    /**
     * @param coverage the coverage to set
     */
    public void setCoverage(int[] coverage) 
    {
        this.coverage = coverage;
    }

    
    /**
     * @return the mutations
     */
    public int[][] getMutations() 
    {
        return mutations;
    }

    
    /**
     * @param mutations the mutations to set
     */
    public void setMutations(int[][] mutations) 
    {
        this.mutations = mutations;
    }

    
    /**
     * @return the insertions
     */
    public int[][] getInsertions() 
    {
        return insertions;
    }

    
    /**
     * @param insertions the insertions to set
     */
    public void setInsertions(int[][] insertions) 
    {
        this.insertions = insertions;
    }

    
    /**
     * @return the deletions
     */
    public int[][] getDeletions() 
    {
        return deletions;
    }

    
    /**
     * @param deletions the deletions to set
     */
    public void setDeletions(int[][] deletions) 
    {
        this.deletions = deletions;
    }

   
    /**
	 * 
	 */
	public void clusterProfile()
	{
	    //Initialize vectors
	    int lengthSequence = sequence.length() != 0 ? sequence.length() : 1;
	    
	    //System.out.println("Cluster Profile: Sequence " + this.sequence + "\tLength " + lengthSequence);
	    
	    this.coverage = new int[lengthSequence];
	    this.mutation = new int[lengthSequence];
	    this.undefined = new int[lengthSequence];
	    this.mutations = new int[12][lengthSequence];
	    this.deletions = new int[4][lengthSequence];
	    this.insertions = new int[4][lengthSequence];
	    
	    
	    //Sequence values
	    String sequence_read;
	    int sequenceIndex;
	    int sequenceOcurrences;
	    String cigar;
	    int length;
	    
	    //Cluster
	    int i = 0; //Initial index
	    int clusterIndex; //
	
	    if(this.sequence.compareTo("") != 0)
	    {
	        for(int h = 0; h < this.sequences.size(); h++) //Move through reads
	        {
	            //Get read information
	            cigar = this.sequences.get(h).getCigar();
	            sequence_read = this.sequences.get(h).getSequence();
	            sequenceOcurrences = this.sequences.get(h).getOccurrences();
	            i = 0; //Initial index
	            
	            //Index
	            clusterIndex = Math.abs(this.sequences.get(h).getStart() - this.getMinPosition());
	            
	            sequenceIndex = 0;
	            
	            
	            for(int k = 0; k < cigar.length(); k++)
	            {
	                if(i == cigar.length() || clusterIndex >= this.sequence.length() || sequenceIndex >= sequence_read.length())
	                	break;
	                
	                if(cigar.charAt(k) == 'M')
	                {
	                	length = Integer.parseInt(cigar.substring(i, k));
	                    
	                    
	                    for(int l = 0; l < length; l++)
	                    {
	                    	if(clusterIndex >= this.sequence.length() || sequenceIndex >= sequence_read.length())
	    	                	break;
	                    	
	                    	this.coverage[clusterIndex] += sequenceOcurrences;
	                    	
	                    	/**
                    	 	 * 0: A - C
                             * 1: A - G
                             * 2: A - T
                             * 3: C - A
                             * 4: C - G
                             * 5: C - T
                             * 6: G - A
                             * 7: G - C
                             * 8: G - T
                             * 9: T - A
                             * 10: T - C
                             * 11: T - G
                            */
	
	                    	if(sequence_read.charAt(sequenceIndex) == 'N')
                            {
	                    		this.undefined[clusterIndex] += sequenceOcurrences;	
	                    	}
	                    	else
	                    	{	
	                    		if(this.sequence.charAt(clusterIndex) == 'A')
	                            {	
	                                if(sequence_read.charAt(sequenceIndex) == 'C')
	                                {
	                                	this.mutations[0][clusterIndex] += sequenceOcurrences;
	                                	this.numberMutations += sequenceOcurrences;
	    	                            this.mutation[clusterIndex] += sequenceOcurrences;
	    	                        }
	                                else
	                                {
	                                    if(sequence_read.charAt(sequenceIndex) == 'G')
	                                    {
	                                    	this.mutations[1][clusterIndex] += sequenceOcurrences;
	                                    	this.numberMutations += sequenceOcurrences;
	        	                            this.mutation[clusterIndex] += sequenceOcurrences;
	        	                        }
	                                    else
	                                    {
	                                        if(sequence_read.charAt(sequenceIndex) == 'T')
	                                        {
	                                        	this.mutations[2][clusterIndex] += sequenceOcurrences;
	                                        	this.numberMutations += sequenceOcurrences;
	            	                            this.mutation[clusterIndex] += sequenceOcurrences;
	            	                        }
	                                    }
	                                }
	                            }
	                            else
	                            {	
	                                if(this.sequence.charAt(clusterIndex) == 'C')
	                                {	
	                                    if(sequence_read.charAt(sequenceIndex) == 'A')
	                                    {
	                                    	this.mutations[3][clusterIndex] += sequenceOcurrences;
	                                    	this.numberMutations += sequenceOcurrences;
	        	                            this.mutation[clusterIndex] += sequenceOcurrences;
	        	                        }
	                                    else
	                                    {	
	                                        if(sequence_read.charAt(sequenceIndex) == 'G')
	                                        {
	                                        	this.mutations[4][clusterIndex] += sequenceOcurrences;
	                                        	this.numberMutations += sequenceOcurrences;
	            	                            this.mutation[clusterIndex] += sequenceOcurrences;
	            	                        }
	                                        else
	                                        {	
	                                            if(sequence_read.charAt(sequenceIndex) == 'T')
	                                            {
	                                            	this.mutations[5][clusterIndex] += sequenceOcurrences;
	                                            	this.numberMutations += sequenceOcurrences;
	                	                            this.mutation[clusterIndex] += sequenceOcurrences;
	                	                        }
	                                        }
	                                    }
	                                }
	                                else
	                                {
	                                    if(this.sequence.charAt(clusterIndex) == 'G')
	                                    {
	                                    	if(sequence_read.charAt(sequenceIndex) == 'A')
	                                        {
	                                    		this.mutations[6][clusterIndex] += sequenceOcurrences;
	                                        	this.numberMutations += sequenceOcurrences;
	            	                            this.mutation[clusterIndex] += sequenceOcurrences;
	            	                        }
	                                        else
	                                        {
	                                            if(sequence_read.charAt(sequenceIndex) == 'C')
	                                            {
	                                            	this.mutations[7][clusterIndex] += sequenceOcurrences;
	                                            	this.numberMutations += sequenceOcurrences;
	                	                            this.mutation[clusterIndex] += sequenceOcurrences;
	                	                        }
	                                            else
	                                            {
	                                                if(sequence_read.charAt(sequenceIndex) == 'T')
	                                                {
	                                                	this.mutations[8][clusterIndex] += sequenceOcurrences;
	                                                	this.numberMutations += sequenceOcurrences;
	                    	                            this.mutation[clusterIndex] += sequenceOcurrences;
	                    	                        }
	                                            }
	                                        }
	                                    }
	                                    else
	                                    {  	
	                                        if(this.sequence.charAt(clusterIndex) == 'T')
	                                        {	
	                                            if(sequence_read.charAt(sequenceIndex) == 'A')            
	                                            {
	                                            	this.mutations[9][clusterIndex] += sequenceOcurrences;
	                                            	this.numberMutations += sequenceOcurrences;
	                	                            this.mutation[clusterIndex] += sequenceOcurrences;
	                	                        }
	                                            else
	                                            { 	
	                                                if(sequence_read.charAt(sequenceIndex) == 'C')
	                                                {
	                                                	this.mutations[10][clusterIndex] += sequenceOcurrences;
	                                                	this.numberMutations += sequenceOcurrences;
	                    	                            this.mutation[clusterIndex] += sequenceOcurrences;
	                    	                        }
	                                                else
	                                                {
	                                                    if(sequence_read.charAt(sequenceIndex) == 'G')
	                                                    {
	                                                    	this.mutations[11][clusterIndex] += sequenceOcurrences;
	                                                    	this.numberMutations += sequenceOcurrences;
	                        	                            this.mutation[clusterIndex] += sequenceOcurrences;
	                        	                        }
	                                                }
	                                        	}
	                                        }
	                                    }
	                                }
	                            }
	                    	}
	                        
	                    	clusterIndex++; //Move in cluster 
	                        sequenceIndex++; //Move in sequence
	                    }
	
	                    i = k + 1; //Move initial index
	                }
	                else
	                {
	                    if(cigar.charAt(k) == 'I')
	                    {
	                        length = Integer.parseInt(cigar.substring(i, k));
	
	                        for(int l = 0; l < length; l++)
	                        {
	                            switch(sequence_read.charAt(sequenceIndex))
	                            {
	                                case 'A': {  this.insertions[0][clusterIndex - 1] += sequenceOcurrences;  }
	                                break;
	                                case 'C': {  this.insertions[1][clusterIndex - 1] += sequenceOcurrences;  }
	                                break;    
	                                case 'G': {  this.insertions[2][clusterIndex - 1] += sequenceOcurrences;  }
	                                break;    
	                                case 'T': {  this.insertions[3][clusterIndex - 1] += sequenceOcurrences;  }
	                                break;
	                            }
	
	                            sequenceIndex++; //Move in sequence
	                        }
	
	                        i = k + 1; //Move initial index
	                    }
	                    else
	                    {
	                        if(cigar.charAt(k) == 'D')
	                        {
	                            length = Integer.parseInt(cigar.substring(i, k));
	                            
	                            for(int l = 0; l < length; l++)
	                            {
	                            	this.coverage[clusterIndex] += sequenceOcurrences;
	                            	
	                            	switch(this.sequence.charAt(clusterIndex))
	                                {
	                                    case 'A': {  this.deletions[0][clusterIndex] += sequenceOcurrences;  }
	                                    break;
	                                    case 'C': {  this.deletions[1][clusterIndex] += sequenceOcurrences;  }
	                                    break;    
	                                    case 'G': {  this.deletions[2][clusterIndex] += sequenceOcurrences;  }
	                                    break;    
	                                    case 'T': {  this.deletions[3][clusterIndex] += sequenceOcurrences;  }
	                                    break;
	                                }
	
	                                clusterIndex++; //Move in cluster 
	                                
	                                if(clusterIndex >= this.sequence.length())
	                                {
	                                	break;
	                                }	
	                            }
	
	                            i = k + 1; //Move initial index
	                        }
	                    }
	                }   
	            }
	        }   
	    } 
	    
	    this.numberMutationsRevelant = this.numberMutations;
	}


	/**
     * 
     * @param position 
     */
    public void compareSNP(int position)
    {
        if(position >= this.minPosition && position < this.maxPosition)
    	{
    		if(this.mutation[position - this.minPosition] != 0)
    		{
    			this.numberMutationsRevelant -= this.mutation[position - this.minPosition];
    			this.mutation[position - this.minPosition] = 0;
    			
    			for(int i = 0; i < 12; i++)
    			{
    				this.mutations[i][position - this.minPosition] = 0;
    			}
    		}
    	}
	}
    
    
    /**
     * 
     * @param id
     * @return
     */
    public String toString(String id)
    {
        String response = "";
        
        response += this.chromosome + "\t" + this.minPosition + "\t" + this.maxPosition + "\t" + id + "\t" + this.getAmountSequences() + "\t" + this.strand;
        response += "\t" + this.sequence + "\t" + this.numberMutations + "\t" + this.numberMutationsRevelant + "\t";

            
        if(sequence.length() != 0)
        {
            for(int i = 0; i < this.coverage.length; i++)
            {    
                response += " " + this.coverage[i];
            }

            response += "\t";
            
            for(int i = 0;i < this.mutation.length; i++)
            {    
                response += " " + this.mutation[i];
            }

            for(int i = 0; i < this.mutations.length; i++)
            {
                response += "\t";
                for(int j = 0; j < this.mutations[0].length; j++)
                    response += " " + this.mutations[i][j];
            }

            for(int i = 0; i < this.insertions.length; i++)
            {
                response += "\t";
                for(int j = 0; j < this.insertions[0].length; j++)
                    response += " " + this.insertions[i][j];
            }

            for(int i = 0; i < this.deletions.length; i++)
            {
                response += "\t";
                for(int j = 0; j < this.deletions[0].length; j++)
                    response += " " + this.deletions[i][j];
            }
            
            response += "\t";
            for(int i = 0; i < this.undefined.length; i++)
            {    
                response += " " + this.undefined[i];
            }
        } 
                
        return response + "\n";
    }
 
    
    /**
     * 
     * @return 
     */
    public String toString_SAM(String id)
    {
        String response = "";
        
        // ExperimentID   Strand   Chromosome  Start  QMap  CIGAR  -  -  -  Sequence Quality  -  -  NumberMutations
        response += id + "\t" + (this.strand.equals("+") ? "0" : "16") + "\t" + this.chromosome + "\t" + this.minPosition + "\t" + this.getClusterQMap();
        response += "\t" + ( sequence.length() != 0 ? sequence.length() : (this.maxPosition - this.minPosition + 1))  + "M" + "\t*\t0\t0";
        response += "\t" + this.sequence + "\t";
        
        for(int i = 0; i < ( sequence.length() != 0 ? sequence.length() : (this.maxPosition - this.minPosition + 1)); i++)
            response += "I";
        
        response += "\tXA:i:_ \tMD:Z:_" + "\tNM:i:" + this.numberMutations + "\n";
        
        return response;
    }
    
    
    /**
     * 
     * @return 
     */
    public String toString_BED(String id)
    {
        String response = "";
        
        //Chromosome    Start   End   Name    Score  Strand
        response += this.chromosome + "\t" + this.minPosition + "\t" + this.maxPosition + "\t" + id + "\t" + this.numberMutations + "\t" + this.strand + "\n";
        
        return response;
    }
    
              
    /**
     * 
     * @return
     */
    public String toString_Fasta()
    {
    	return ">" + this.chromosome + ":" + this.minPosition + ":" + this.maxPosition + ":" + this.strand + "\n" + this.sequence + "\n";
    }
 
    
    
    public String toString_Prior(int window)
    {
    	double prior[] = new double[this.mutations[0].length];
    	String result = "";
    	double total = 0.0;
    	
    	for(int i = 0; i < prior.length; i++)
    	{
    		for(int j = 0; j < window; j++)
    		{
    			prior[i] += (i + j) < prior.length ? (double)this.mutation[i] : 0.0;
    		}
    		
    		total += prior[i];
    	}
    	
    	result += ">" + this.chromosome + ":" + this.minPosition + ":" + this.maxPosition + ":" + this.strand + "\n";
    	
    	for(int i = 0; i < prior.length; i++)
    	{
    		prior[i] /= total;
    		result += (Double.isNaN(prior[i]) ? "0.0" : prior[i]) + " ";
    	}
    	
    	return  result + "\n";
    }
   
    
    /**
     * 
     * @return 
     */
    public int getClusterQMap()
    {
        int average = 0;
        
        for(int i = 0; i < this.sequences.size(); i++)
            average += this.sequences.get(i).getMap_quality();
        
        
        average /= this.sequences.size();
        
        return average;
    }
    
    
    /**
     * 
     * @return 
     */
    public String toStringReads()
    {
        String response = "";
        
        for(int i = 0; i < this.sequences.size(); i++)
        {
            response += this.sequences.get(i).toString() + "\n";
        }
        
        return response;
    }

    
    /**
     * 
     * @return 
     */
    public String toStringSequences()
    {
        String response = "\n\n\n";
        int spaces;
        
        response += this.sequence + "\n";
        
        for(int i = 0; i < this.sequences.size(); i++)
        {
            spaces = this.sequences.get(i).getStart() - this.minPosition;
            
            for(int j = 0; j < spaces; j++)
            {
                response += " ";
            }
            
            response += this.sequences.get(i).getSequence() + "\n";
        }
        
        return response;
    }
    
    
    public String getLog()
    {
    	return this.logSNP;
    }
    
    
    public boolean SNPcontact()
    {
    	return this.numberMutations == this.numberMutationsRevelant ? true : false;
    }
    
    
    public void getClusterSequence(String fastaPath)
    {
         String line = "";
         String sequence = "";
         int index = 0;
         String path = fastaPath + this.chromosome + ".fa";
         
         try 
         {
             File fasFile = new File (path);
             
             if(fasFile.exists())
             {
                 BufferedReader br = new BufferedReader(new FileReader(fasFile));
                 line = br.readLine(); //Header
                 line = br.readLine(); //First line
                 
                 boolean successful = false;
                 int start, end;
                 
                 start = this.getMinPosition();
 				 end = this.getMaxPosition();
                 
 				 sequence = "";
 				 successful = false;
 					
 				 do
 	            {
 	                if((index + 50) > (start - 1))
 	                {
 	                    if((index + 50) > end)
 	                    {
 	                    	if(!(((start - 1) % 50) >= (end % 50)))
 	                    	{
 	                    		sequence += line.substring((start - 1) % 50, end % 50);
 	                    
 	                    	}
 	                    }
 	                    else
 	                    {
 	                    	sequence += line.substring((start - 1) % 50, 50);
 	                    	index += 50;
 	                    	line = br.readLine(); //Next line
 	                    	
 	                    	while((index + 50) <= end)
 	                        {
 	                    		sequence += line;
 	                    		line = br.readLine(); //Next line
 	                    		index += 50;
 	                        }
 	                           
 	                    	sequence += line.substring(0, (end % 50));
 	                    }
 	                    successful = true;
 	                }
 	                
 	                if(!successful)
 	                {
 	                	index += 50;
 	                	line = br.readLine(); //Next line
 	                }
 	            }
 	            while(line != null && !successful);
 		                
 				 if(successful)
 					this.setSequence(sequence.toUpperCase()); 
 				 
 				 br.close();
             }	 
         } 
         catch (FileNotFoundException ex)   {} 
         catch (IOException ex)   {}
     }
    
    
    /**
     * 
     * @param fastaPath
     */
    public void sequenceValidation(String fastaPath)
    {
    	if(sequence.length() != (this.maxClusterLength - this.minClusterLength + 1))
    	{
    		this.getClusterSequence(fastaPath);
    	}
    }
    
 
    /**
     * 
     * @return
     */
    public int getMutationsRanking()
    {
    	int counter = 0;
    	
    	for(int i = 0; i < this.mutations[10].length; i++)
    	{
    		counter += this.mutations[10][i];
    	}
    	
    	return counter;
    }
    
    /**
     * 
     * @return
     */
    public int getDensityRanking()
    {
    	int counterCoverage = 0;
    	int counterMutations = 0;
    	
    	
    	for(int i = 0; i < this.coverage.length; i++)
    	{
    		counterCoverage += this.coverage[i];
    	}
    	
    	for(int i = 0; i < this.mutations[10].length; i++)
    	{
    		counterMutations += this.mutations[10][i];
    	}
    	
    	return counterMutations / counterCoverage;
    }
    
    
    public void change(Cluster cluster)
    {
    	//Attributes
        this.sequences = cluster.getSequences(); //Sequences inside the cluster
        this.minPosition = cluster.getMinPosition(); // Start position of the cluster
        this.maxPosition = cluster.getMaxPosition(); //End position of the cluster
        this.strand = cluster.getStrand(); //Strand associated at the cluster
        this.chromosome = cluster.getChromosome(); //Chromosome associated at the cluster
    	
        //General Indicators
        this.minClusterLength = cluster.getMinClusterLength(); //Minimum length of a cluster
        this.maxClusterLength = cluster.getMaxClusterLength(); //Maximum length of a cluster
        this.minClusterSequences = cluster.getMinClusterSequences(); //Minimum number of sequences in a cluster
        this.maxClusterSequences = cluster.getMaxClusterSequences(); //Maximum number of sequences in a cluster
        this.numberMutations = 0; //NUmber of Mutations
        this.numberMutationsRevelant = 0; //Number of mutations with SNPs filter
        
        //Profile
        this.sequence = cluster.getSequence();  //Cluster sequence
        this.coverage = cluster.getCoverage(); //Coverage vector
        //this.mutation = cluster.getMutation(); //General mutations vector
        //this.undefined = cluster.getUndefined(); //General mutations vector
        this.mutations = cluster.getMutations(); //Mutations vector by mutation (A-C, A-G, A-T, ..., T-G)
        this.insertions = cluster.getInsertions(); //Insertions vector by nucleotide base (A, C, G, T)
        this.deletions = cluster.getDeletions(); //Deletions vector by nucleotide base (A, C, G, T)
    }
}