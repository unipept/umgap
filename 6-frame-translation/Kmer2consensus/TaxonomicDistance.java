	
public class TaxonomicDistance{
	public static void main(String[] args){
		String trueLineage = args[0];
		String lineage = args[1];
		String[] trueLin = trueLineage.split(";");
		String[] lin = lineage.split(";");
		int agreement = 0;
		while(agreement < lin.length && agreement < trueLin.length && lin[agreement].trim().equals(trueLin[agreement].trim())){	
            		agreement +=1;
        	}
		int disagreement = lin.length + trueLin.length - 2*agreement;
		System.out.println(disagreement);
	}
}
