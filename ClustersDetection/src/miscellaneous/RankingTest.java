package miscellaneous;

import java.util.Vector;

import structures.Cluster;

public class RankingTest {
	
	public Vector<Cluster> test = new Vector<Cluster>();
	
	public RankingTest() {}


	public void insert(Cluster cluster)
	{
		int mutations = cluster.getMaxClusterLength();
		int amountSequences = cluster.getAmountSequences(); 
		boolean success = false;
		boolean notFound = false;
		
		int middle;
		int min = 0;
		int max = this.test.size() - 1;
		
		while(!notFound && !success)
		{
			middle = (max + min) / 2;
			
			if(mutations == test.get(middle).getMaxClusterLength())
			{
				//TODO order by amountSequences
				
				if(amountSequences == test.get(middle).getAmountSequences())
				{
					//TODO order by ClusterLength
				}
				else
				{
					if(amountSequences < test.get(middle).getAmountSequences())
					{
						do
						{
							middle--;
						}
						while(amountSequences < test.get(middle).getAmountSequences());
						
					}
					else
					{
						do
						{
							middle++;
						}
						while(amountSequences > test.get(middle).getAmountSequences());
					}	
					
				}
			}
			else
			{
				if(mutations < test.get(middle).getMaxClusterLength())
				{
					max = middle;
				}
				else
				{
					min = middle;
				}
			}	
		}
	}
}
