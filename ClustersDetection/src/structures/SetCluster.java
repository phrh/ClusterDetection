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
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Vector;
import java.util.logging.Level;
import java.util.logging.Logger;

import miscellaneous.Functions;
import filesProcess.CD_FromSam;


/**
 *
 * @author Eng. Paula Reyes, Ph.D. -- M.Sc. Eng. Carlos Sierra
 * Universidad Antonio Narino 
 */
public class SetCluster extends Thread {
    
    //Attributes
    private Vector<Cluster> set = null; //Clusters without filters

    private int minClusterLength; //Minimum length of a cluster
    private int maxClusterLength; //Maximum length of a cluster
    private int minClusterSequences; //Minimum number of sequences in a cluster
    private int maxClusterSequences; //Maximum number of sequences in a cluster
    private String chromosome; //Chromosome associated at setCluster
    
    private String folderResults; //Folder to save the results of the set of clusters
    private String prefixFile; //
    private int minLength; //
    private int minSequences; //
    private boolean saveDetails; //
    private String fastaPath; //
    private String snpPath; //
    
    
    private String log = "";
    
    
    /**
     * Zero-parameter constructor
     */
    public SetCluster(){}
  
    
    /**
     * 
     * @param folderResults
     * @param prefixFile
     * @param minLength
     * @param minSequences
     * @param saveDetails
     * @param chromosome
     * @param fastaPath
     * @param snpPath
     * @param startPositionPath
     * @param startMutationPath
     * @param otherPositionPath
     * @param otherMutationPath
     */
    public SetCluster(String folderResults, String prefixFile, int minLength, int minSequences, boolean saveDetails, String chromosome, String fastaPath, String snpPath, String startPositionPath, String startMutationPath, String otherPositionPath, String otherMutationPath)
    {
        this.set = new Vector<Cluster>();
        this.minClusterLength = 1000000;
        this.minClusterSequences = 1000000;
        this.maxClusterLength = 0;
        this.maxClusterSequences = 0;
        
        this.folderResults = folderResults;
        this.prefixFile = prefixFile;
        this.minLength = minLength;
        this.minSequences = minSequences;
        this.saveDetails = saveDetails;
        this.chromosome = chromosome;
        this.fastaPath = fastaPath;
        this.snpPath = snpPath;
        //this.startPositionPath = startPositionPath;
        //this.startMutationPath = startMutationPath;
        //this.otherPositionPath = otherPositionPath;
        //this.otherMutationPath = otherMutationPath;
    }

    
    /**
     * @return the set
     */
    public Vector<Cluster> getSet() 
    {
        return set;
    }

    /**
     * @param set the set to set
     */
    public void setSet(Vector<Cluster> set) 
    {
        this.set = set;
    }
    
    
    
    @Override
    public void run()
    {
    	long startTime = System.nanoTime(); //Start time counter
    	double elapsedTimeInSec;
    	
    	File path = new File(this.folderResults);
        if (!path.exists())  // if folder doesn't exists, then create it
            path.mkdirs();
        
        //log += "Clean overlapping...\n";
        this.cleanOverlapping();
        
        elapsedTimeInSec = (System.nanoTime() - startTime) * 1.0e-9;
        log += "Elapsed Time: " + elapsedTimeInSec + " seconds\n";
        
        if(this.set.size() > 0)
        {
        	log += "\nProfile clusters...\n";
        	this.profiles();
            elapsedTimeInSec = (System.nanoTime() - startTime) * 1.0e-9;
            log += "Elapsed Time: " + elapsedTimeInSec + " seconds\n";

            log += "\nProcess SNPs...\n";
            String temp_bed = this.chromosome.substring(0, 3) + "_" + this.chromosome.substring(3);
        	this.processSNP(this.snpPath + "bed_" + temp_bed + ".bed");
        	elapsedTimeInSec = (System.nanoTime() - startTime) * 1.0e-9;
            log += "Elapsed Time: " + elapsedTimeInSec + " seconds\n";
        	
            log +=  "\nSaving clusters...\n";
            this.saveClusters();  // TODO   
        }
        else
        {
        	log +=  "No clusters to process...\n\n";
        }
        
        elapsedTimeInSec = (System.nanoTime() - startTime) * 1.0e-9;
        log += "Elapsed Total Time: " + elapsedTimeInSec + " seconds\n";
        
        interrupt();
    }

    
    /**
     * 
     * @param startPosition
     * @param endPosition
     * @param strand
     * @param chromosome
     * @return
     */
    public int getClusterIndex(int startPosition, int endPosition, String strand, String chromosome) 
    {
        int min = 0; //Current lower limit of search segment
    	int max = this.set.size(); //Current upper limit of search segment
    	int size = this.set.size(); //Size of the cluster
 		int middle; //Split point of search segment
 		int index = -1;
		
 		
 		do //Use binary search to find the respective index
        {
            middle = (min + max) / 2;

            if(middle >= 0 && (middle < size) && this.set.get(middle).isSequence(startPosition, endPosition, strand, chromosome))
            {
                index = middle; //The cluster exists
                break;
            }
            else
            {
            	if (middle >= 0 && middle < size && startPosition > this.set.get(middle).getMinPosition())
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

       /**
        * This method take reported SnPs and compare with clusters
        * @param snpPath
        */
	   private void processSNP(String snpPath)
	   {
	       int i = 0; //Clusters index
	       int positionStart;
	       int positionEnd;
	       String strand;
	       String chromosome;
	       String line = "";
	       String[] content;
	       
	        try 
	        {
	            File bedFile = new File (snpPath);
	            
	            if(bedFile.exists())
	            {
	            	/*File fileSNP = new File(this.startPositionPath);
		        	FileWriter fwSNP;
		            if (!fileSNP.exists()) 
		            {
		                try 
		                {
		                	fileSNP.createNewFile();
		                } 
		                catch (IOException ex) 
		                {
		                    Logger.getLogger(Chromosome.class.getName()).log(Level.SEVERE, null, ex);
		                    System.out.println("\n\nFile " + fileSNP);
		                    System.exit(0);
		                }
		            } 
		            

		            File fileoSNP = new File(this.otherPositionPath);
		        	FileWriter fwoSNP;
		            if (!fileSNP.exists()) 
		            {
		                try 
		                {
		                	fileoSNP.createNewFile();
		                } 
		                catch (IOException ex) 
		                {
		                    Logger.getLogger(Chromosome.class.getName()).log(Level.SEVERE, null, ex);
		                    System.out.println("\n\nFile " + fileoSNP);
		                    System.exit(0);
		                }
		            } 

		            
		            File fileSNPdis = new File(this.startMutationPath);
		        	if (!fileSNPdis.exists()) 
		            {
		                try 
		                {
		                	fileSNPdis.createNewFile();
		                } 
		                catch (IOException ex) 
		                {
		                    Logger.getLogger(Chromosome.class.getName()).log(Level.SEVERE, null, ex);
		                    System.out.println("\n\nFile " + fileSNPdis);
		                    System.exit(0);
		                }
		            }
		        	
		        	File fileoSNPdis = new File(this.otherMutationPath);
		        	if (!fileoSNPdis.exists()) 
		            {
		                try 
		                {
		                	fileoSNPdis.createNewFile();
		                } 
		                catch (IOException ex) 
		                {
		                    Logger.getLogger(Chromosome.class.getName()).log(Level.SEVERE, null, ex);
		                    System.out.println("\n\nFile " + fileoSNPdis);
		                    System.exit(0);
		                }
		            }
		        	*/
		            
	                BufferedReader br = new BufferedReader(new FileReader(bedFile));
	                line = br.readLine(); //Header
	                line = br.readLine(); //First line
	
	                //fwSNP = new FileWriter(fileSNP.getAbsoluteFile(), true);
	    			//BufferedWriter bwSNP = new BufferedWriter(fwSNP);
	    			
	    			//fwoSNP = new FileWriter(fileoSNP.getAbsoluteFile(), true);
	    			//BufferedWriter bwoSNP = new BufferedWriter(fwoSNP);
	                
	    			while(line != null)
	                {
	                    content = line.split("\t");
	                    positionStart = Integer.parseInt(content[1]);
	                    positionEnd = Integer.parseInt(content[2]);
	                    
	                    if((positionEnd - positionStart) == 1)
	                    {
	                    	strand = content[4];
	                    	chromosome = content[0];
	                    	
	                    	i = this.getClusterIndex(positionEnd, positionEnd, strand, chromosome);
	                    	
		                    if(i >= 0 && i < this.set.size())
		                    {
		                        this.set.get(i).compareSNP(positionEnd);
		                        
		                        /*if(this.set.get(i).getMinPosition() == positionEnd)
		                        {
		                        	bwSNP.write(this.set.get(i).toString_BED("0"));
		                        }
		                        else
		                        {
		                        	bwoSNP.write(this.set.get(i).toString_BED("1"));  //TODO 
		                        }*/
		                    }
	                    }
		                    
	                    line = br.readLine(); //Next line
	                }
	                
	    			//bwSNP.close();
	    			//bwoSNP.close();
	                br.close();
	            }    
	        } 
	        catch (FileNotFoundException ex) 
	        {
	            Logger.getLogger(SetCluster.class.getName()).log(Level.SEVERE, null, ex);
	        } 
	        catch (IOException ex) 
	        {
	            Logger.getLogger(SetCluster.class.getName()).log(Level.SEVERE, null, ex);
	        }
	    }


/**
    * This method takes the clusters set given, and clean the overlapping clusters merge them. 
    * Additionally, here we make a statistical process about the number of sequences and length (maximum minus minimum position) of each cluster.
    * That information is used to build histograms and another statistical tools.
    * @param clustersSet
    * @param prefixeFile
    * @return clusters-set
    */
   private void cleanOverlapping()
   {
       int minPosition, maxPosition; //Current minimum and maximum position of the respective cluster
       String strand; //Current strand of the respective cluster
       String chromosome; //Current chromosome of the respective cluster
       int negativeStrand = 0; //Number of clusters with negative strand
       int positiveStrand = 0; //Number of clusters with positive strand
       boolean overlaps = false;
       
       //int temp = 0;
       
       //Ask in all the current clusters
       do
       {
	       overlaps = false;
    	   for(int i = this.set.size() - 1; i >= 0; i--) //Move up to down - Last to First
	       {
	           minPosition = this.set.get(i).getMinPosition(); //Get minimum position of the cluster
	           maxPosition = this.set.get(i).getMaxPosition(); //Get maximum position of the cluster
	           strand = this.set.get(i).getStrand();  //Get strand of the cluster
	           chromosome = this.set.get(i).getChromosome();
	
	           for(int j = 1; j <= 5; j++)
	           {
	        	   if((i - j) >= 0 && this.set.get(i - j).isSequence(minPosition, maxPosition, strand, chromosome)) //Verify overlapping with the adjacent cluster
	        	   {
	        		   this.mergeClusters(i - j, i); //Is overlapping, then merge clusters
	        		   overlaps = true;
	        		   break;
	        	   }
	           }
	       }
	   }
       while(overlaps);
       
       //Obtain statistics
       //System.out.println("Size:" + this.set.size() + "\t temp " + temp);
       for(int i = this.set.size() - 1; i >= 0; i--) 
       {
    	   if(this.set.get(i).getAmountSequences() < 2)
    	   {
    		   this.set.remove(i);
    	   }
    	   else
    	   {
    		    //Get boundaries
	            this.minClusterLength = (this.set.get(i).length() < this.minClusterLength) ? this.set.get(i).length() : this.minClusterLength;
	            this.maxClusterLength = (this.set.get(i).length() > this.maxClusterLength) ? this.set.get(i).length() : this.maxClusterLength;
	            this.minClusterSequences = (this.set.get(i).getAmountSequences() < this.minClusterSequences) ? this.set.get(i).getAmountSequences() : this.minClusterSequences;
	            this.maxClusterSequences = (this.set.get(i).getAmountSequences() > this.maxClusterSequences) ? this.set.get(i).getAmountSequences() : this.maxClusterSequences;
	       }
       }
       
       int[] occurrences = new int[this.maxClusterSequences + 1]; 
       int[] lenghts = new int[this.maxClusterLength + 1]; 
       
       //Obtain statistics
       for(int i = 0; i < this.set.size(); i++) 
       {
    	   occurrences[this.set.get(i).getAmountSequences()] += 1;
           lenghts[this.set.get(i).length()] += 1;

           if(this.set.get(i).getStrand().equals("+"))
        	   positiveStrand++;
           else
        	   negativeStrand++;
       }
       
       log += "Clusters with Positive Strand: " + positiveStrand + "\tClusters with Negative Strand: " + negativeStrand + "\n";

       
       //---------- OCURRENCES ----------//
        String fileOcurrences = this.folderResults + this.prefixFile + "Ocurrences.txt";
        File fileSaveO = new File(fileOcurrences);
        
        FileWriter fwO;
        if (!fileSaveO.exists()) 
        {
            try 
            {
                fileSaveO.createNewFile();
            } 
            catch (IOException ex) 
            {
                Logger.getLogger(Chromosome.class.getName()).log(Level.SEVERE, null, ex);
                System.out.println("\n\nFile " + fileOcurrences);
                System.exit(0);
            }
        }

       try 
       {
           fwO = new FileWriter(fileSaveO.getAbsoluteFile(), true);
           BufferedWriter bwO = new BufferedWriter(fwO);
           
           for(int j = 0; j < occurrences.length; j++)
           {
               if(occurrences[j] != 0) //If there's some occurrence of j-number of sequences
               {
                   //Clusters    AmountSequences
                   bwO.write(j + " \t " + occurrences[j] + "\n");
               }
           }
           
           bwO.close();
       } 
       catch (IOException ex) 
       {
           Logger.getLogger(CD_FromSam.class.getName()).log(Level.SEVERE, null, ex);
           System.exit(0);
       }
        
        
        
        //---------- LENGTHS ----------//
        String fileLengths = this.folderResults + this.prefixFile + "Lengths.txt";
        File fileSaveL = new File(fileLengths);
        
        FileWriter fwL;
        if (!fileSaveL.exists()) 
        {
            try 
            {
                fileSaveL.createNewFile();
            } 
            catch (IOException ex) 
            {
                Logger.getLogger(Chromosome.class.getName()).log(Level.SEVERE, null, ex);
                System.out.println("\n\nFile " + fileLengths);
                System.exit(0);
            }
        }

       try 
       {
           fwL = new FileWriter(fileSaveL.getAbsoluteFile(), true);
           BufferedWriter bwL = new BufferedWriter(fwL);
           
           for(int j = 0; j < lenghts.length; j++)
           {
               if(lenghts[j] != 0) //If there's some occurrence of j-length of cluster
               {
                   //Clusters    AmountSequences
                   bwL.write(j + " \t " + lenghts[j] + "\n");
               }
           }
           
           bwL.close();
       } 
       catch (IOException ex) 
       {
           Logger.getLogger(CD_FromSam.class.getName()).log(Level.SEVERE, null, ex);
           System.exit(0);
       }
    }


   /**
    * This method merge the clusters in the given indexes (Cluster-indexToMerge inside Cluster-indexMain).
    * @param indexMain
    * @param indexToMerge
    * @return clusters set
    */
   private void mergeClusters(int indexMain, int indexToMerge)
   {
       //Extract one per one the clusters of "indexToMerge"
       for(int i = 0; i < this.set.get(indexToMerge).getSequences().size(); i++)
       {
           this.set.get(indexMain).addSequence(this.set.get(indexToMerge).getSequence(i));
       }

       //End to copy all the sequences
       this.set.remove(indexToMerge);
    } 


   /**
    * 
    */
   public void getClusterSequence()
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
                
                for(int i = 0; i < this.set.size(); i ++)
				{
					start = this.set.get(i).getMinPosition();
					end = this.set.get(i).getMaxPosition();
                
					sequence = "";
					successful = false;
					
					if(this.set.get(i).getStrand().equals("+"))
	                {
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
							this.set.get(i).setSequence(sequence.toUpperCase()); 
	                } 
				}
                
                br.close();
            
                
                br = new BufferedReader(new FileReader(fasFile));
                line = br.readLine(); //Header
                line = br.readLine(); //Read first line
        		
        		index = 0;
                
                for(int i = 0; i < this.set.size(); i ++)
				{
					start = this.set.get(i).getMinPosition();
					end = this.set.get(i).getMaxPosition();
                
					sequence = "";
					successful = false;
					
					if(this.set.get(i).getStrand().equals("-"))
	                {
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
							this.set.get(i).setSequence(sequence.toUpperCase()); 
	                }
				}
                
                br.close();
                
                //TODO Validation
            }
        } 
        catch (FileNotFoundException ex)   {} 
        catch (IOException ex)   {}
    }
   
   
   /**
    * 
    */
   private void profiles()
   {
       this.getClusterSequence();
	   
	   for(int i = 0; i < this.set.size(); i++) //Move up to down - Last to First
       {
            this.set.get(i).clusterProfile();
       }
   }
   
   
   /**
    * 
     * @param clusters
     */
    private void saveClusters()
    {
       String file = this.folderResults + this.prefixFile + "Clusters_" + this.minLength + ".0.txt";
       String fileFilterAmountSequences = this.folderResults +  this.prefixFile + "Clusters_" + this.minLength + "." + this.minSequences + ".txt";
       String path = this.folderResults + this.prefixFile + "Details/";

       File folderPath = new File(path);
       if (!folderPath.exists() && this.saveDetails) 
           folderPath.mkdirs();

       File fileClusters = new File(file);
       FileWriter fwClusters;
       if (!fileClusters.exists()) 
       {
           try 
           {
               fileClusters.createNewFile();
           } 
           catch (IOException ex) 
           {
               Logger.getLogger(Chromosome.class.getName()).log(Level.SEVERE, null, ex);
               System.out.println("\n\nFile " + file);
               System.exit(0);
           }
       }
       
       File fileFilter = new File(fileFilterAmountSequences);
       FileWriter fwFilter;
       if (!fileClusters.exists()) 
       {
           try 
           {
               fileFilter.createNewFile();
           } 
           catch (IOException ex) 
           {
               Logger.getLogger(Chromosome.class.getName()).log(Level.SEVERE, null, ex);
               System.out.println("\n\nFile " + fileFilterAmountSequences);
               System.exit(0);
           }
       }
       
       
       int contTotal = 0; //Number of clusters
       int contFiltered = 0; //Number of clusters after amount sequences filter
       
       try 
       {
           
    	   
    	   fwClusters = new FileWriter(fileClusters.getAbsoluteFile(), true);
           BufferedWriter bwClusters = new BufferedWriter(fwClusters);
           
           fwFilter = new FileWriter(fileFilter.getAbsoluteFile(), true);
           BufferedWriter bwFilter = new BufferedWriter(fwFilter);
       
           
	       for(int i = 0; i < this.set.size(); i++) 
	       {
	    	   contTotal++; 
		
	           	//ID    Chromosome    Start   End   Name    Score  Strand
	            bwClusters.write(contTotal + "\t" + this.set.get(i).toString(this.chromosome + "." + contTotal));
	            //this.log += this.set.get(i).getLog();
	
	            if(this.set.get(i).getAmountSequences() >= this.minSequences)
	            {
	                contFiltered++; 
	
	                //ID   Chromosome    Start   End   Name    Score  Strand
	                bwFilter.write(contFiltered + "\t" + this.set.get(i).toString(this.chromosome + "." + contFiltered));
	                    
	                if(this.saveDetails)
	                	saveEachCluster(this.set.get(i), contFiltered, this.folderResults + this.prefixFile + "Details/Cluster");
	            }
	       }
	       
	       bwClusters.close();
	       bwFilter.close();
	       
	       for(int i = this.set.size() - 1; i >= 0; i--) 
	       {
	    	   if(this.set.get(i).getAmountSequences() < this.minSequences)
	           {
	        	  this.set.remove(i); 
	           }
	       } 
    
	       
	   } 
       catch (IOException ex) 
       {
           Logger.getLogger(CD_FromSam.class.getName()).log(Level.SEVERE, null, ex);
           System.exit(0);
       }
       
       
       log += "Clusters: " + contTotal + "\tClusters Filtered by MinSequences: " + contFiltered + "\n";
   }


   /**
    * This method save in a file (plain text) the information of the cluster given.
    * @param cluster
    * @param id
    * @param prefixFile 
    */
   private void saveEachCluster(Cluster cluster, int id, String prefixFile)
   {
	   //Name of file to save cluster
       String file = prefixFile + "_" + id + "__" + cluster.getMinPosition() + "-" + cluster.getMaxPosition() + ".txt";
       File fileSave = new File(file);
       FileWriter fwSaveClusters;
       if (!fileSave.exists()) 
       {
           try 
           {
               fileSave.createNewFile();
           } 
           catch (IOException ex) 
           {
               Logger.getLogger(Chromosome.class.getName()).log(Level.SEVERE, null, ex);
               System.out.println("\n\nFile " + file);
               System.exit(0);
           }
       }

       try 
       {
    	   	fwSaveClusters = new FileWriter(fileSave.getAbsoluteFile(), true);
    	   	BufferedWriter bwSaveClusters = new BufferedWriter(fwSaveClusters);
    	       
    	   		
    	   	//Chromosome    Start   End   Name   Occurrences(Score)  Strand
    	   	for(int i = 0; i < cluster.getSequences().size(); i++)
    	   	{
    	   		bwSaveClusters.write(cluster.getSequences().get(i).toString() + "\n");
    	   	}
               
            bwSaveClusters.write("\n\n\n" + cluster.getSequence() + "\n");
            
            int spaces;
            String response;
            
            for(int i = 0; i < cluster.getSequences().size(); i++)
            {
               spaces = cluster.getSequences().get(i).getStart() - cluster.getMinPosition();
               response = ""; 
               
               for(int j = 0; j < spaces; j++)
               {
                   response += " ";
               }
               
               response += cluster.getSequences().get(i).getSequence() + "\n";
               bwSaveClusters.write(response);
           }
               
           bwSaveClusters.close();
       } 
       catch (IOException ex) 
       {
          Logger.getLogger(CD_FromSam.class.getName()).log(Level.SEVERE, null, ex);
          System.exit(0);
       }
   }
   
   /**
    * This method saves information in Chromosome's log.
    */
   public void saveLog()
   {
       Functions.saveInFile(this.folderResults + "/Log.txt", this.log);
   }
}