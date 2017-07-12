/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package drawlca4fgs;

import java.io.File;
import java.io.FileReader;
import java.io.StringReader;
import java.util.List;
import java.util.Scanner;
import com.opencsv.CSVReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Methode analoog aan ScoreReads, maar dan om een extra stukje tikz-code te genereren 
 * dat een plot maakt van een opsplitsing van een FGR-voorspelling
 * @author Aranka
 */
public class DrawLca4FGS {
    private static int depth;
    private static double unit;
    private static double begin;
    /**
     * @param args the command line arguments
     * 0: hight at which to draw
     * 1: start position
     * 2: total length
     * 3: k = 0 for tryptic peptides k > 0 for k-mer of length k
     * 4: the protein
     * 5: name of / path to the file with the lca's
     * 6: name of / path to the file with the lineages of the lca
     * 7: true lineage
     */
    public static void main(String[] args) {
        depth = Integer.parseInt(args[0]);
        begin = Double.parseDouble(args[1]);
        unit = (double) 24 / (double) Integer.parseInt(args[2]);
        int k = 0;
        boolean tp;
        if (Integer.parseInt(args[3]) == 0 ){
            tp = true;
        }else{
            tp = false;
            k = Integer.parseInt(args[3]);
        }
        String protein = args[4];
        System.out.println("\\draw (0,"+ (depth+1) + ") -- (24,"+ (depth+1) + ") node[right] {FGS+};");
        System.out.println("\\draw[RoyalBlue, line width = 4pt, line cap = round] ("+((double)begin*unit)+","+(depth+1)+") -- ("+((double)(begin+protein.length())*unit )+","+(depth+1)+");");
        try {
            Scanner lca = new Scanner(new File(args[5]));
            String lcaHeader = lca.nextLine();
            CSVReader lineageR = new CSVReader(new StringReader(""));
            String[] trueLineage = new String[0];
            if(args.length > 6){
                File lineageF = new File(args[6]);
                lineageR = new CSVReader(new FileReader(lineageF),';');
                trueLineage = args[7].split(";");
            }
            if (tp){
                ArrayList<Peptide> peptides = new ArrayList<>();
                while(lca.hasNextLine()){
                    String nextLCA = lca.nextLine();
                    if(args.length > 6){
                        if(nextLCA.contains("root")){
                            String[] root = new String[1];
                            root[0] = "root";
                            peptides.add(new Peptide(nextLCA, root));
                        }else{
                            String[] lineage = lineageR.readNext();
                            peptides.add(new Peptide(nextLCA, lineage));;
                        }
                    }else{
                        peptides.add(new Peptide(nextLCA));
                    }        
                }
//                printProteins((double)24/(double)frame.length(),proteins);
                if(args.length > 6){
                    drawKnownFrame(protein,depth,peptides,true,trueLineage);
                }
                drawUnknownFrame(protein,depth+0.5,peptides,true);
            }else{
    //      Specifiek voor K-meren
                ArrayList<Kmer> kmers = new ArrayList<>();
                while(lca.hasNextLine()){
                    String nextLCA = lca.nextLine();
                    if(args.length > 6){
                        if(nextLCA.contains("root")){
                            String[] root = new String[1];
                            root[0] = "root";
                            kmers.add(new Kmer(nextLCA,k, root));;
                        }else{
                            String[] lineage = lineageR.readNext();
                            kmers.add(new Kmer(nextLCA,k, lineage));;
                        }
                    }else{
                        kmers.add(new Kmer(nextLCA,k));
                    }
                }
                    int[] frameDepths = getCoverageDepth(protein,kmers,k);
                    if(args.length > 6){
                        plotCoverageDepths(1,depth+2,frameDepths,k,2);
                        drawUnknownKmerFrame(protein,depth+1,true,kmers,k);
                        drawKnownKmerFrame(protein, depth, true,kmers, k, trueLineage);
                    }else{
                        plotCoverageDepths(1,depth+1,frameDepths,k,4);
                        drawUnknownKmerFrame(protein, depth, true ,kmers,k);
                    }
                }
        } catch (FileNotFoundException ex) {
            Logger.getLogger(DrawLca4FGS.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(DrawLca4FGS.class.getName()).log(Level.SEVERE, null, ex);
        }
        
    }
    public static void drawUnknownFrame(String frame, double framenr, List<Peptide> peptides, boolean forward){
        for(Peptide pept:peptides){
            double start;
            double end;
            if(forward){
                start = (double) (frame.indexOf(pept.aminoSeq)+begin)*unit;
                end = start + (double) pept.length*unit;
            }else{
                end = (double)frame.length()*unit-(double)frame.indexOf(pept.aminoSeq)*unit-(double)begin*unit;
                start = end - (double) pept.length*unit;    
            }
            String taxonName = pept.taxonName;
            if(taxonName.equalsIgnoreCase("root")){
                System.out.println("\\draw[root] ("+start+","+framenr+") -- ("+end+","+framenr+");");
            }else{
                System.out.println("\\draw["+pept.taxonRank+"] ("+start+","+framenr+") -- ("+end+","+framenr+");");
            }
        }
    }
    
    public static void drawKnownFrame(String frame, int framenr, List<Peptide> peptides, boolean forward, String[] trueLineage){
        String lineageString="";
        for(int i=0; i<trueLineage.length;i++){
            lineageString=lineageString + trueLineage[i];
        }
        for (Peptide pept:peptides){
            double start;
            double end;
            if(forward){
                start = (double) (begin+frame.indexOf(pept.aminoSeq))*unit;
                end = start + (double) pept.length*unit;
            }else{
                end = (double)frame.length()*unit-(double)frame.indexOf(pept.aminoSeq)*unit-(double)begin*unit;
                start = end - (double) pept.length*unit;    
            }
            String taxonName = pept.taxonName;
            if(taxonName.equalsIgnoreCase("root")){
                System.out.println("\\draw[root] ("+start+","+framenr+") -- ("+end+","+framenr+");");
            }else{
                if(lineageString.contains(taxonName)){
                    System.out.println("\\draw["+pept.taxonRank+"] ("+start+","+framenr+") -- ("+end+","+framenr+");");
                }else{
                    String[] lineage = pept.lineage;
                    printDisagreement(lineage,trueLineage,start,end,framenr);
                }
            }
        }
    }
    
    public static void printDisagreement(String[] lineage, String[] trueLin, double start, double end, int framenr){
        int[] disagreement = new int[2];
        int agreement = 0;
        while(agreement< lineage.length && agreement< trueLin.length &&lineage[agreement].trim().equals(trueLin[agreement])){
            agreement +=1;
        }
        disagreement[0] = lineage.length + 1 + trueLin.length - 2*agreement;
        disagreement[1] = 100-disagreement[0];
        System.out.println("\\draw[red!"+disagreement[1]+", line width = 4pt, line cap = round] ("+start+","+framenr+") -- ("+end+","+framenr+");");
        System.out.println("\\node[below, font=\\small] at ("+(start+end)/(double)2+","+framenr+") {"+disagreement[0]+"};");
    }
    
    public static void printKmerDisagreement(String[] lineage, String[] trueLin, double start, int framenr){
        int[] disagreement = new int[2];
        int agreement = 0;
        while(agreement< lineage.length && agreement< trueLin.length &&lineage[agreement].trim().equals(trueLin[agreement])){
            agreement +=1;
        }
        disagreement[0] = lineage.length + 1 + trueLin.length - 2*agreement;
        disagreement[1] = 100-disagreement[0];
        double above = ((double)framenr + 0.2)/2;
        double below = ((double)framenr - 0.2)/2;
        System.out.println("\\draw[red!"+disagreement[1]+", line width = 2pt, line cap = butt] ("+start+","+above+") -- ("+start+","+below+");");
        System.out.println("\\node[draw=none, below, font=\\sffamily\\fontsize{2}{1}\\selectfont] at ("+start+","+(below+0.2)+") {"+disagreement[0]+"};");
    }
    
    private static int[] getCoverageDepth(String frame, List<Kmer> kmers, int k){
        int[] depths = new int[frame.length()];
        int cur_pos = 0;
        for(Kmer kmer:kmers){
            int start = (int) begin + frame.indexOf(kmer.aminoSeq, cur_pos);
            cur_pos = start;
            for(int i=start;i<start+k;i++){
                depths[i] = depths[i] + 1;
            }
        }
        return depths;
    }
    
    public static void plotCoverageDepths(int frame, int diepte, int[] depths, int k, int factor){
        int n = depths.length;
        for(int i=0;i<depths.length;i++){
            int fact = depths[i]*100/k;
            double above = ((double)diepte + 0.2)/factor;
            double below = ((double)diepte - 0.2)/factor;
            if(frame<4){
                System.out.println("\\draw [RoyalBlue!"+fact+",line width=2pt] ("+(i)*unit+","+above+") -- ("+(i)*unit+","+below+");");
            }else{
                System.out.println("\\draw [RoyalBlue!"+fact+",line width=2pt] ("+(n-i)*unit+","+above+") -- ("+(n-i)*unit+","+below+");");
            }
        }
    }
    public static void drawUnknownKmerFrame(String frame, int framenr, boolean forward, List<Kmer> kmers, int k){
        int pos=0;
        for(Kmer kmer: kmers){
            double start;
            start = begin + frame.indexOf(kmer.aminoSeq, pos);
            pos = (int) start + 1;
            start = (double) start*unit;
            if(! forward){
                start = (double)frame.length()*unit - start;   
            }
            String taxonName = kmer.taxonName;
            double above = ((double)framenr + 0.2)/2;
            double below = ((double)framenr - 0.2)/2;
            if(taxonName.equalsIgnoreCase("root")){
                System.out.println("\\draw[root, line cap=butt, line width = 2pt] ("+start+","+above+") -- ("+start+","+below+");");
            }else{
                System.out.println("\\draw["+kmer.taxonRank+", line cap=butt, line width = 2pt] ("+start+","+above+") -- ("+start+","+below+");");
            }
        }
    }
    public static void drawKnownKmerFrame(String frame, int framenr, boolean forward, List<Kmer> kmers, int k, String[] trueLineage){
        String lineageString="";
        for(int i=0; i<trueLineage.length;i++){
            lineageString=lineageString + trueLineage[i];
        }
        int pos=0;
        for(Kmer kmer: kmers){
            double start;
            start = begin + frame.indexOf(kmer.aminoSeq, pos);
            pos = (int) start + 1;
            start = (double) start*unit;
            if(! forward){
                start = (double)frame.length()*unit - start;   
            }
            String taxonName = kmer.taxonName;
            double above = ((double)framenr + 0.2)/2;
            double below = ((double)framenr - 0.2)/2;
            if(taxonName.equalsIgnoreCase("root")){
                System.out.println("\\draw[root, line cap=butt, line width = 2pt] ("+start+","+above+") -- ("+start+","+below+");");
            }else{
                if(lineageString.contains(taxonName)){
                    System.out.println("\\draw["+kmer.taxonRank+", line cap=butt, line width = 2pt] ("+start+","+above+") -- ("+start+","+below+");");
                }else{
                    String[] lineage = kmer.lineage;
                    printKmerDisagreement(lineage,trueLineage,start,framenr);
                }
            }
        }
    }
    
}
