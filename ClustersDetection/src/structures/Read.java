package structures;

/**
 * This class represents the behavior of a sequence.
 * @author Eng. Paula Reyes, Ph.D. -- M.Sc. Eng. Carlos Sierra
 * Antonio Narino University
 */
public class Read {

    //Attributes
    private String firstData; //Experiment ID
    private String chromosome; //Chromosome name
    private String strand; //1 if strand equals to +; 0 if strand equals to -
    private int map_quality; //255 indicates that isn't available quality
    private int start; //Start position of the sequence
    private int end;  //Start position of the sequence
    private String cigar; //CIGAR of sequence: Matchs, Deletions, and Insertions
    private double score; //Number of mutations
    private int occurrences; //Occurrences of the sequence inside the cluster
    private String sequence; //RNA fragment
    private String quality; //Quality string


    /**
     * Zero-parameters constructor. Default values for attributes. 
     */
    public Read() 
    {
        this.firstData = "";
        this.chromosome = "";
        this.strand = "";
        this.start = 0;
        this.end = 0;
        this.map_quality = -1;
        this.score = 0;
        this.sequence = "";
        this.quality = "";
        this.cigar = "";
        this.occurrences = 0;
    }


    /**
     * Constructor with parameters (Recommended)
     * @param experiment_id
     * @param chromosome
     * @param strand
     * @param start
     * @param end
     * @param qMap
     * @param score
     * @param sequence
     * @param quality 
     * @param cigar
     */
    public Read(String firstData, String chromosome, String strand, int start, int end, int qMap, double score, String sequence, String quality, String cigar) 
    {
        this.firstData = firstData;
        this.chromosome = chromosome;
        this.strand = strand;
        this.start = start;
        this.end = end;
        this.map_quality = qMap;
        this.score = score;
        this.sequence = sequence;
        this.quality = quality;
        this.cigar = cigar;
        this.occurrences = 1;
    }


    /**
     * This method returns the experiment-id of the sequence.
     * @return experiment-id
     */
    public String getFirstData() 
    {
        return firstData;
    }


    /**
     * This method sets the experiment-id of the sequence.
     * @param id of firstData to set
     */
    public void setFirstData(String firstData) 
    {
        this.firstData = firstData;
    }


    /**
     * This method returns the sequence start position.
     * @return start position
     */
    public int getStart() 
    {
        return this.start;
    }


    /**
     * This method sets the start position value of the sequence.
     * @param start the start position to set
     */
    public void setStart(int start) 
    {
        this.start = start;
    }


    /**
     * This method returns the sequence end position.
     * @return end position
     */
    public int getEnd() 
    {
        return this.end;
    }


    /**
     * This method sets the end position value of the sequence.
     * @param end the end position to set
     */
    public void setEnd(int end) 
    {
        this.end = end;
    }

    
    /**
     * @return the cigar
     */
    public String getCigar() 
    {
        return cigar;
    }

    
    /**
     * @param cigar the cigar to set
     */
    public void setCigar(String cigar) 
    {
        this.cigar = cigar;
    }


    /**
     * This method returns the chromosome associated at the sequence.
     * @return the chromosome
     */
    public String getChromosome() 
    {
        return this.chromosome;
    }


    /**
     * This method sets the chromosome associated at the sequence.
     * @param chromosome the chromosome to set
     */
    public void setChromosome(String chromosome) 
    {
        this.chromosome = chromosome;
    }


    /**
     * This method return the strand of the sequence.
     * @return strand
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
     * @return the score
     */
    public double getScore() 
    {
        return this.score;
    }


    /**
     * @param score the score to set
     */
    public void setScore(double score) 
    {
        this.score = score;
    }


    /**
     * This method returns the occurrences of the read in the cluster.
     * @return occurrences
     */
    public int getOccurrences() 
    {
        return this.occurrences;
    }


    /**
     * This method sets the occurrences of the read in a cluster.
     * @param occurrences the occurrences to set
     */
    public void setOccurrences(int occurrences) 
    {
        this.occurrences = occurrences;
    }

    
    /**
     * This method returns the mapping quality.
     * @return the map_quality
     */
    public int getMap_quality() 
    {
        return map_quality;
    }

    
    /**
     * This method sets the mapping quality.
     * @param map_quality the map_quality to set
     */
    public void setMap_quality(int map_quality) 
    {
        this.map_quality = map_quality;
    }

    
    /**
     * This method returns the read sequence.
     * @return the sequence
     */
    public String getSequence() 
    {
        return sequence;
    }

    
    /**
     * This method sets the sequence read.
     * @param sequence the sequence to set
     */
    public void setSequence(String sequence) 
    {
        this.sequence = sequence;
    }

    
    /**
     * This method returns the sequence quality.
     * @return sequence quality
     */
    public String getQuality() 
    {
        return quality;
    }

    
    /**
     * This method sets the sequence quality.
     * @param quality the quality to set
     */
    public void setQuality(String quality) 
    {
        this.quality = quality;
    }
	
	
    /**
     * Increase the occurrences of sequence inside the cluster.
     */
    public void increaseOccurrence()
    {
        this.occurrences++;
    }
    
    
    @Override
    public String toString()
    {
        String response = "";
        
        response += this.firstData + "\t" + this.chromosome + "\t" +  this.start + "\t" + this.end + "\t NameX \t" + this.occurrences + " \t" + this.cigar + " \t" + this.strand + " \t" + this.sequence + " \t" + this.quality + " i\t" + this.score;
        
        return response;
    }
}