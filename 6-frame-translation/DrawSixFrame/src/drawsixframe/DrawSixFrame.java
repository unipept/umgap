/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package drawsixframe;

import com.opencsv.CSVReader;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.TreeSet;


/**
 *
 * @author Aranka
 */
public class DrawSixFrame {
    private static final String correctTaxonomy = "cellular organisms; Bacteria; Proteobacteria; Gammaproteobacteria; Pseudomonadales; Moraxellaceae; Acinetobacter; Acinetobacter calcoaceticus/baumannii complex; Acinetobacter baumannii";
    private static final TreeSet<String> taxonomy = new TreeSet<>();
    private static final double pagewidth = (double) 24;
    private static final ArrayList<String> proteinSeq = new ArrayList<>();
    
    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        proteinSeq.add("MKISDLMTYHGCKNRKELSEKTGYSTVTLWKWENNGIPARTQAVLQVKTKGKLKADLQALTA");
        proteinSeq.add("MSLHSRIRQKLEEKKLRAADLARATKKSPVAVKKWLDGTSVPTAENLKVIAKFLGVSDDWLLYGGPVEQESNNLPQLNVLDIEAFKQKYNIPDSEDAVKFVQTSDKPFPIQKRYVPVKAYSKMGMDGYFTDMGYEGNAGDGYVPTHSAGPRAYGIKGTGDSMFPAIRNGWYVVCDPDAELVPTEFVQVCLKDGRCTIKEFIGINNDVLSLIAVNGGERLSFNMDEVESITAITDIVPPSQHRQEHPYSH");
        taxonomy.addAll(Arrays.asList(correctTaxonomy.split("; ")));
        CSVReader reader;
        try {
            reader = new CSVReader(new FileReader("found_lcas_double.txt"));
            BufferedReader sixframe = new BufferedReader(new FileReader("double.sixframe"));
            List<String[]> lcas = reader.readAll();
            Iterator<String[]> it = lcas.iterator();
            it.next();
            String header;
            int framenr = 1;
            String[] nextline = it.next();
            String frame = "";
            double unit = 1;
            while((header=sixframe.readLine())!=null){
                framenr+=1;
                frame = sixframe.readLine();
                unit = pagewidth/(double)frame.length();
                for(String seq: proteinSeq){
                    if(frame.contains(seq)){
                        double proteinstart;
                        double proteinstop;
                        if(framenr>=2 && framenr <5){
                            proteinstart = (double) frame.indexOf(seq)*unit;
                            proteinstop = proteinstart + (double) seq.length()*unit;
                        }else{
                            proteinstop = (double)frame.length()*unit-(double)frame.indexOf(seq)*unit;
                            proteinstart = proteinstop - (double) seq.length()*unit;
                        }
                    System.out.println("\\draw[protein] ("+proteinstart+",0) -- ("+proteinstop+",0);");
                    }
                }
                boolean sameHeader = nextline[0].equals(header);
                while(it.hasNext() && sameHeader){
                    if(framenr>=2 && framenr <5){
                        double start = (double) frame.indexOf(nextline[1])*unit;
                        double end = start + (double) nextline[1].length()*unit;
                        String taxonName = nextline[3];
                        if(taxonName.equalsIgnoreCase("root")){
                            System.out.println("\\draw[root] ("+start+",-"+framenr+") -- ("+end+",-"+framenr+");");
                        }else{
                            if(taxonomy.contains(taxonName)){
                                System.out.println("\\draw["+nextline[4]+"] ("+start+",-"+framenr+") -- ("+end+",-"+framenr+");");
                            }else{
                                System.out.println("\\draw[wrong] ("+start+",-"+framenr+") -- ("+end+",-"+framenr+");");
                            }
                        }
                    }else{
                        double end = (double)frame.length()*unit-(double)frame.indexOf(nextline[1])*unit;
                        double start = end - (double) nextline[1].length()*unit;
                        String taxonName = nextline[3];
                        if(taxonName.equalsIgnoreCase("root")){
                            System.out.println("\\draw[root] ("+start+",-"+framenr+") -- ("+end+",-"+framenr+");");
                        }else{
                            if(taxonomy.contains(taxonName)){
                                System.out.println("\\draw["+nextline[4]+"] ("+start+",-"+framenr+") -- ("+end+",-"+framenr+");");
                            }else{
                                System.out.println("\\draw[wrong] ("+start+",-"+framenr+") -- ("+end+",-"+framenr+");");
                            }
                        }
                    }
                    nextline = it.next();
                    sameHeader = nextline[0].equals(header);
                }
            }
            double end = (double)frame.length()*unit-(double)frame.indexOf(nextline[1])*unit;
            double start = end - (double) nextline[1].length()*unit;
            String taxonName = nextline[3];
            if(taxonName.equalsIgnoreCase("root")){
                System.out.println("\\draw[root] ("+start+",-"+framenr+") -- ("+end+",-"+framenr+");");
            }else{
                if(taxonomy.contains(taxonName)){
                    System.out.println("\\draw["+nextline[4]+"] ("+start+",-"+framenr+") -- ("+end+",-"+framenr+");");
                }else{
                    System.out.println("\\draw[wrong] ("+start+",-"+framenr+") -- ("+end+",-"+framenr+");");
                }
            }    
            } catch (IOException ex) {
            System.out.println("File found_lcas.txt not found");
        }
    }
    
}
