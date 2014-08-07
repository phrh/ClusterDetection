/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package miscellaneous;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author Eng. Paula Reyes, Ph.D. -- M.Sc. Eng. Carlos Sierra
 * Universidad Antonio Narino 
 */
public class ChromosomeSubsequence 
{
    public ChromosomeSubsequence()
    {
        
    }
    
    
    public void getSequence()
    {
        int start = 103471450;
        int end =  103471500;
        String sequence = "";
        String chrm = "chr1";
        String line = "";
        int index = 0;
        System.out.println("Length: " + (end - start + 1));
        
        try 
        {
            BufferedReader br = new BufferedReader(new FileReader(chrm + ".fa"));
            line = br.readLine(); //Header
            line = br.readLine(); //Next line
            boolean successful = false;
            
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
				sequence = sequence.toUpperCase(); 
            
            br.close();
            System.out.println("Seq: " + sequence + "\tLength:" + sequence.length());
            
        } 
        catch (FileNotFoundException ex) 
        {
            Logger.getLogger(ChromosomeSubsequence.class.getName()).log(Level.SEVERE, null, ex);
        } 
        catch (IOException ex) 
        {
            Logger.getLogger(ChromosomeSubsequence.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    

    public static void main(String args[])
    {
        ChromosomeSubsequence cs = new ChromosomeSubsequence();
        cs.getSequence();
    }
}
