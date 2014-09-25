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

import java.io.File;
import java.util.Vector;
import java.util.logging.Level;
import java.util.logging.Logger;

import miscellaneous.Functions;


/**
 * This class represents a general structure of a chromosome by SAM, when just use a little of information of the chromosome
 * @author Eng. Paula Reyes, Ph.D. -- M.Sc. Eng. Carlos Sierra
 * Universidad Antonio Narino 
 */
public class Chromosome extends Thread
{
    //Attributes
    private String name; //Name of the chromosome
    private int chrLenght; //Length of the chromosome
    private SetCluster clusters = null; //Clusters without filters
    private SetCluster clustersFilter = null; //Clusters with read-length filter
    
    //Parameters for process
    private String folderResults;
    private int minLength; 
    
    /**
     * Zero-parameters constructor of the class
     */
    public Chromosome()  
    {
        //Define attributes
        this.name = "unname";
        this.chrLenght = 0;
        this.clusters = new SetCluster(); 
        this.clustersFilter = new SetCluster(); 
        
        //Define parameters
        this.folderResults = "";
        this.minLength = 0;
    } 
    
    
    /**
     * Constructor with parameters of the class
     */
    public Chromosome(String name, int chrLenght, String folderResults, int minSequences, int minLength, String fastaPath, String snpPath)
    {
        //Define attributes
        this.name = name;
        this.chrLenght = chrLenght;
        this.clusters = new SetCluster(folderResults + "/" + this.name, "/Total_", 0, 0, false, this.name, fastaPath, snpPath, folderResults + "/SNP_startPosition.txt", folderResults + "/SNP_affectMutations.txt", folderResults + "/SNP_otherPosition.txt", folderResults + "/SNP_otherAffectMutations.txt"); 
        this.clustersFilter = new SetCluster(folderResults + "/" + this.name, "/Filter_", minLength, minSequences, true, this.name, fastaPath, snpPath, folderResults + "/SNP_filter_startPosition.txt", folderResults + "/SNP_filter_affectMutations.txt", folderResults + "/SNP_filter_otherPosition.txt", folderResults + "/SNP_filter_otherAffectMutations.txt"); 
        
        //Define parameters
        this.folderResults = folderResults + "/" + this.name;
        this.minLength = minLength;
        
        File folderPath = new File(this.folderResults);
        if (!folderPath.exists()) 
            folderPath.mkdirs();
    }
    
    
    /**
     * 
     * @return 
     */
    public String getChromosomeName()
    {
        return this.name;
    }
    
    
    /**
     * 
     * @param name 
     */
    public void setChromosomeName(String name)
    {
        this.name = name;
    }

        
    /**
     * 
     * @return 
     */
    public int getChrLength()
    {
        return this.chrLenght;
    }
    
    
    /**
     * 
     * @param chrLenght 
     */
    public void setChrLength(int chrLenght)
    {
        this.chrLenght = chrLenght;
    }

    
    /**
     * @return the clusters
     */
    public SetCluster getClusters() 
    {
    	return clusters;
    }

    
    /**
     * @param clusters the clusters to set
     */
    public void setClusters(SetCluster clusters) 
    {
        this.clusters = clusters;
    }

    
    /**
     * @return the clustersFilter
     */
    public SetCluster getClustersFilter() 
    {
        return clustersFilter;
    }

    
    /**
     * @param clustersFilter the clustersFilter to set
     */
    public void setClustersFilter(SetCluster clustersFilter) 
    {
        this.clustersFilter = clustersFilter;
    }
    
    
    @Override
    public void run()
    {
        this.clusters.start();
        this.clustersFilter.start();
        
        try 
        {
            Thread.sleep(100);
        } 
        catch (InterruptedException ex) 
        {
            Logger.getLogger(Chromosome.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        try 
        {
            Functions.saveInFile(folderResults + "/Log.txt", "\n\nTOTAL:\n");
            this.clusters.join();
            this.clusters.saveLog();
            
            Functions.saveInFile(folderResults + "/Log.txt", "\n\nFILTERED:\n");
            this.clustersFilter.join();
            this.clustersFilter.saveLog();
            
            Functions.saveInFile(folderResults + "/Log.txt", "\nFinish chromosome");
            
            System.out.println("\nReady  " + this.name);
            this.interrupt();
        } 
        catch (InterruptedException ex) 
        {
            Logger.getLogger(Chromosome.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    
 
    /**
     * 
     * @param id_experiment
     * @param chromosome
     * @param start
     * @param end
     * @param strand
     * @param qmap
     * @param numberMutations
     * @param sequence
     * @param quality
     * @param minLength 
     */
    public void addSequence(String id_experiment, String chromosome, int start, int end, String strand, int qmap, int numberMutations, String sequence, String quality, String cigar)
    {
        int indexCluster; //Index where insert sequence
        int indexClusterFilter; //Index where insert sequence in filter vector (If the sequence has no contraints problems)
        
        //Search possible cluster
        indexCluster = this.insertSequence(start, end, strand, this.clusters.getSet(), chromosome);

        //Validate if the cluster not exists
        if(indexCluster == this.clusters.getSet().size() || !this.clusters.getSet().get(indexCluster).isSequence(start, end, strand, chromosome)) 
            this.clusters.getSet().add(indexCluster, new Cluster(start, end, strand, chromosome)); //Create new cluster
        
        //Add sequence to cluster without filters
        this.clusters.getSet().get(indexCluster).addSequence(new Read(id_experiment, chromosome, strand, start, end, qmap, numberMutations, sequence, quality, cigar));
       
        
        //Validate if the sequence pass the filter -> Minimum length minLength bases
        if(Math.abs(end - start) > this.minLength)
        {
            //Search possible cluster
            indexClusterFilter = this.insertSequence(start, end, strand, this.clustersFilter.getSet(), chromosome);

            //Validate if the cluster not exists
            if(indexClusterFilter == this.clustersFilter.getSet().size() || !this.clustersFilter.getSet().get(indexClusterFilter).isSequence(start, end, strand, chromosome)) 
                this.clustersFilter.getSet().add(indexClusterFilter, new Cluster(start, end, strand, chromosome)); //New cluster
 
            //Add sequence to cluster without filters
            this.clustersFilter.getSet().get(indexClusterFilter).addSequence(new Read(id_experiment, chromosome, strand, start, end, qmap, numberMutations, sequence, quality, cigar));
        }
    }
   
     
    /**
     * This method return the cluster where the sequence must to be classified.
     * @param startPosition
     * @param endPosition
     * @param strand
     * @param clustersSet
     * @param chromosome
     * @return 
     */
    public int insertSequence(int startPosition, int endPosition, String strand, Vector<Cluster> clustersSet, String chromosome) 
    {
        int min = 0; //Current lower limit of search segment
    	int max = clustersSet.size(); //Current upper limit of search segment
    	int size = clustersSet.size(); //Size of the cluster
 		int middle; //Split point of search segment
 		int index = -1;
			
        do //Use binary search to find the respective index
        {
            middle = (min + max) / 2;

            if(middle >= 0 && (middle < size) && clustersSet.get(middle).isSequence(startPosition, endPosition, strand, chromosome))
            {
                index = middle; //The cluster exists
                break;
            }
            else
            {
                if (middle >= 0 && middle < size && startPosition > clustersSet.get(middle).getMinPosition())
                    min = middle + 1; //Discard the lower part
                else 
                    max = middle; //Discard the lower part
            }
	   		    
            if(min == max)
                index = min;
        }
        while (min < max);
			
        return index;
    }
}