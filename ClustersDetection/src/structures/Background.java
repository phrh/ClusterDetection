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

import java.util.Vector;

/**
 * 
 * @author Eng. Paula Reyes, Ph.D. -- M.Sc. Eng. Carlos Sierra
 * Universidad Antonio Narino
 */
public class Background {

	//Attributes
	private String name;
	private Vector<Read> readsPositive = new Vector<Read>();
	private Vector<Read> readsNegative = new Vector<Read>();
	
	/**
	 * 
	 * @param name
	 */
	public Background(String name) 
	{
		this.name = name;
	}
	
	
	/**
	 * @return the name
	 */
	public String getName() {
		return name;
	}


	/**
	 * @param name the name to set
	 */
	public void setName(String name) {
		this.name = name;
	}


	/**
	 * @return the readsPositive
	 */
	public Vector<Read> getReadsPositive() 
	{
		return readsPositive;
	}


	/**
	 * @param readsPositive the readsPositive to set
	 */
	public void setReadsPositive(Vector<Read> readsPositive) 
	{
		this.readsPositive = readsPositive;
	}


	/**
	 * @return the readsNegative
	 */
	public Vector<Read> getReadsNegative() 
	{
		return readsNegative;
	}


	/**
	 * @param readsNegative the readsNegative to set
	 */
	public void setReadsNegative(Vector<Read> readsNegative) 
	{
		this.readsNegative = readsNegative;
	}


	/**
	 * 
	 * @param start
	 * @param end
	 * @param strand
	 * @param score
	 */
	public void insert(int start, int end, String strand, int score)
	{
		if(strand.equals("+"))
		{
			this.readsPositive.add(new Read("", name, strand, start, end, 0, score, "", "", ""));
		}
		else
		{
			this.readsNegative.add(new Read("", name, strand, start, end, 0, score, "", "", ""));
		}
	}

	
	/**
	 * 
	 * @param start
	 * @param end
	 * @param strand
	 * @return
	 */
	public boolean search(int start, int end, String strand)
	{
		boolean finded = false;
		
		if(strand.equals("+"))
		{
			for(int i = 0; i < readsPositive.size(); i++)
			{
				if(readsPositive.get(i).getStrand().equals(strand))
				{
					if((start <= readsPositive.get(i).getStart() && end >= readsPositive.get(i).getEnd()) || (readsPositive.get(i).getStart() <= start && readsPositive.get(i).getEnd() >= end))
					{
						return true;
					}
				}
			}
		}
		else
		{
			for(int i = 0; i < readsNegative.size(); i++)
			{
				if(readsNegative.get(i).getStrand().equals(strand))
				{
					if((start <= readsNegative.get(i).getStart() && end >= readsNegative.get(i).getEnd()) && (readsNegative.get(i).getStart() <= start && readsNegative.get(i).getEnd() >= end))
					{
						return true;
					}
				}
			}
		}
		
		return finded;
	}
	
}