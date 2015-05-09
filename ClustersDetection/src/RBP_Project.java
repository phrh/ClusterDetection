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


import filesProcess.CD_FromSam;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;

import miscellaneous.Functions;


/**
 *
 * @author Eng. Paula Reyes, Ph.D. -- M.Sc. Eng. Carlos Sierra
 * Universidad Antonio Narino 
 */
public class RBP_Project {
    public static void main(String[] args)
    {  
      	if(args.length > 0)
        {    
            //Parameters of clusterDetections
      		String dataset = "";
      		String experiment = "";
            int minSequences = 5;
            int minLenghts = 20;
            String folder = "results";  
            String mail = "";
            String fastaPath = "";
            String snpPath = "";
            String noisePath = "";
            String[] mutations = null;
            String[] deletions = null;
            String[] insertions = null;
            boolean exportSAM = false;
            boolean exportBED = false;
            boolean filterEntropy = true;
            int window = 8;
            int amount = 1000;
            int bg_score = 3;
            
            
      		System.gc(); // Call Garbage Collector
            
      		String setupFile = args[0]; //File to extract reads
            String line;
            
            try 
			{
            	BufferedReader br;
            	br = new BufferedReader(new FileReader(setupFile));
            	String[] field_value;
            	
            	line = br.readLine();
                
                while(line != null)
                {
                	field_value = line.split("=");
                	
                	
                	//Read parameters
                	switch(field_value[0])
                	{
                		case "DataSet":
                		{
                			String[] datasetPath = field_value[1].split("/");
                			experiment = datasetPath[datasetPath.length - 1];		
                			dataset = field_value[1];
                		}
                		break;
                		case "Filter_Minimum_Ocurrences":
                		{
                			minSequences = Integer.parseInt(field_value[1]);
                		}
                		break;
                		case "Filter_Minimum_Read_length":
                		{
                			minLenghts  = Integer.parseInt(field_value[1]);
                		}
                		break;
                		case "Folder_Results":
                		{
                			 folder = field_value[1];
                		}
                		break;
                		case "E-mail":
                		{
                			mail = field_value[1];
                		}
                		break;
                		case "FASTA_Folder":
                		{
                			fastaPath = field_value[1];
                			
                			File fastaPathF = new File(fastaPath);
                            if (!fastaPathF.exists())
                            {
                            	System.out.println("Genome folder not exists.");
                            	System.exit(0);
                            }
                        }
                		break;
                		case "SNP_Folder":
                		{
                			snpPath = field_value[1];
                			File snpPathF = new File(snpPath);
                            if (!snpPathF.exists())
                            {
                            	System.out.println("SNPs folder not exists.");
                            	System.exit(0);
                            }
                		}
                		break;
                		case "Noise_Reads":
                		{
                			noisePath = field_value[1];
                		}
                		break;
                		case "Export_SAM":
                		{
                			 exportSAM = field_value[1].toLowerCase().equals("yes") ? true : false;
                		}
                		break;
                		case "Export_BED":
                		{
                			exportBED = field_value[1].toLowerCase().equals("yes") ? true : false;
                		}
                		break;
                		case "Mutations":
                		{
                			mutations = field_value[1].split(",");
                		}
                		break;
                		case "Insertions":
                		{
                			insertions = field_value[1].split(",");
                		}
                		break;
                		case "Deletions":
                		{
                			deletions = field_value[1].split(",");
                		}
                		break;
                		case "Filter_Entrophy":
                		{
                			filterEntropy = field_value[1].toLowerCase().equals("yes") ? true : false;
                		}
                		break;
                		case "Window":
                		{
                			window = Integer.parseInt(field_value[1]);
                		}
                		break;
                		case "Ranking":
                		{
                			amount = Integer.parseInt(field_value[1]);
                		}
                		break;
                		case "Background_Score":
                		{
                			bg_score = Integer.parseInt(field_value[1]);
                		}
                		break;
                	}
                	
                	line = br.readLine();
                }
                
                br.close();
			} 
			catch (FileNotFoundException e1) 
			{
				e1.printStackTrace();
			} 
            catch (IOException e) 
			{
				e.printStackTrace();
			}
            
            
            if(dataset.equals(""))
            {
            	System.out.println("DataSet parameter is empty or invalid.");
            	System.exit(0);
            }
            
            //Create folder to save the experiments results
            File folderPath = new File(folder);
            if (!folderPath.exists()) 
                folderPath.mkdirs();
            
            
            //Create the general log file
            String fileLog = folder + "/Log.txt";
           
            
            //Read and process the file
            try 
            {
               CD_FromSam cdSam = new CD_FromSam();
               cdSam.readFile(dataset, folder, minLenghts, minSequences, fileLog, fastaPath, snpPath, filterEntropy, exportSAM, exportBED, noisePath, mutations, deletions, insertions, window, amount, setupFile, bg_score);
                    
               Thread.sleep(500);
               Functions.sentMail(mail, "The experiment is complete.\n", "DataSet: " + experiment);
            } 
            catch (InterruptedException e) 
            {
                Logger.getLogger(RBP_Project.class.getName()).log(Level.SEVERE, null, e);
                Functions.sentMail(mail, "ERROR: Read and process file fail.\n" + e.getMessage() + "\n", "DataSet: " + experiment);
                System.exit(0);
            } 
            catch (FileNotFoundException e) 
            {
                Logger.getLogger(RBP_Project.class.getName()).log(Level.SEVERE, null, e);
                Functions.sentMail(mail, "ERROR: File not found.\n" + e.getMessage() + "\n", "DataSet: " + dataset);
                System.exit(0);
            }
        }
        else
        {
            System.out.println("Incorrect parameters...");
        }
        
        System.gc(); // Call Garbage Collector
     }  
}  