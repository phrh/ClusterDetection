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


package miscellaneous;

import filesProcess.CD_FromSam;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.Properties;
import java.util.logging.Level;
import java.util.logging.Logger;

import javax.mail.Message;
import javax.mail.MessagingException;
import javax.mail.Session;
import javax.mail.Transport;
import javax.mail.internet.InternetAddress;
import javax.mail.internet.MimeMessage;

import structures.Chromosome;


/**
 *
 * @author Eng. Paula Reyes, Ph.D. -- M.Sc. Eng. Carlos Sierra
 * Universidad Antonio Narino 
 */
public class Functions {

    /**
     * 
     * @param file
     * @param text 
     */
    public static void saveInFile(String file, String text)
    {
        File fileSave = new File(file);
        
        //300 segs aprox
        FileWriter fwl;
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
           fwl = new FileWriter(fileSave.getAbsoluteFile(), true);
           try(BufferedWriter bwl = new BufferedWriter(fwl))
           {
                bwl.write(text);
           }
       } 
       catch (IOException ex) 
       {
           Logger.getLogger(CD_FromSam.class.getName()).log(Level.SEVERE, null, ex);
           System.exit(0);
       }
    }
    
    
     /**
     * 
     * @param to
     * @param content
     * @param experimentName 
     */
    public static void sentMail(String to, String content, String experimentName)
    {
        String from = "phreyes@gmail.com"; // Sender E-mail
        String host = "localhost"; //SMTP Host
        Properties properties = System.getProperties(); //System Properties
        properties.setProperty("mail.smtp.host", host); // Setup Mail Server

        // Get Default Session object.
        Session session = Session.getDefaultInstance(properties);

         try
         {
            // Create a default MimeMessage object.
            MimeMessage mail = new MimeMessage(session);

            mail.setFrom(new InternetAddress(from)); //Set "From"
            mail.addRecipient(Message.RecipientType.TO, new InternetAddress(to)); //Set "To"
            mail.setSubject("End Run: " + experimentName); //Set "Subject"

            String timeStamp = new SimpleDateFormat("yyyyMMdd_HHmmss").format(Calendar.getInstance().getTime());
            
            mail.setContent("<h2>" + experimentName + "</h2><br><br>" + content + "<br><br>Date: " + timeStamp, "text/html" );

            Transport.send(mail); //Send message
            System.out.println("Sent mail successfully....");
        }
        catch (MessagingException e) 
        {
            System.out.println("Fail Sent mail.\n" + e.getMessage());
        }
    }
}