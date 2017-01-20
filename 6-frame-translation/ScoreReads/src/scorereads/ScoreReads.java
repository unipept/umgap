/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package scorereads;

import com.opencsv.CSVReader;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Scanner;
import java.util.TreeMap;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author Aranka
 */
public class ScoreReads {
    private static Map<String,Double> taxonomy_score = new TreeMap<>();
    private static String title = "Acinetobacter Baumannii";
    private static String subtitle = "";

    /**
     * @param args the command line arguments
     *              0: =0 for tryptic peptides
     *                 =k for k-mers (k is an integer)
     *              1: Title for the plot
     *              2: name of / path to the file with six frame translation
     *              3: name of / path to the file with the lca's
     *              4: the positions of present proteins
     *              5(optional): name of / path to a file with lineages belonging to the lca's
     *              6(optional): String with the lineage of the organism to which the DNA belongs
     */
    public static void main(String[] args) {
        if (args.length < 5 || args.length > 7){
            System.out.println("Arguments: tp/kmer-indicator    plot-title  sixframe-file   lca-file    proteinpositions");
        }else{
            boolean tp;
            int k = 0;
            if (Integer.parseInt(args[0]) == 0 ){
                tp = true;
            }else{
                tp = false;
                k = Integer.parseInt(args[0]);
            }
            title = args[1];
            printHeader();
            try {
                fillTaxScore("taxonomy_score.txt"); // zorg gewoon dat de gewenste score in een bestand met deze naam zit
            } catch (FileNotFoundException ex) {
                Logger.getLogger(ScoreReads.class.getName()).log(Level.SEVERE, null, ex);
            }
            File sixFrame = new File(args[2]);
            String[] proteins = new String[0];
            if (! args[4].equals(":")){
                proteins = args[4].split(":");
            }
            try {
                BufferedReader sixframe = new BufferedReader(new FileReader(sixFrame));
                Scanner lca = new Scanner(new File(args[3]));
                CSVReader lineageR = new CSVReader(new StringReader(""));
                String[] trueLineage = new String[0];
                if(args.length > 5){
                    File lineageF = new File(args[5]);
                    lineageR = new CSVReader(new FileReader(lineageF),';');
                    trueLineage = args[6].split(";");
                }
                if(tp){
                    lca.nextLine();
                    printEmptyLines(2);
                    if (args.length > 5){
                        printEmptyLines(8);
                    }
                }
                String lcaHeader = lca.nextLine();
                int diepte = 2;
                int forward = 1;
                String header;
                while((header=sixframe.readLine())!=null){
                    String frame = sixframe.readLine();
    //      Specifiek voor Triptische Peptiden  
                    if (tp){
                        ArrayList<Peptide> peptides = new ArrayList<>();
                        while(lca.hasNextLine() && lcaHeader.contains(header)){
                            if(args.length > 5){
                                if(lcaHeader.contains("root")){
                                    String[] root = new String[1];
                                    root[0] = "root";
                                    peptides.add(new Peptide(lcaHeader, taxonomy_score, root));
                                }else{
                                    String[] lineage = lineageR.readNext();
                                    peptides.add(new Peptide(lcaHeader, taxonomy_score, lineage));
                                }
                            }else{
                                peptides.add(new Peptide(lcaHeader,taxonomy_score));
                            }
                            lcaHeader = lca.nextLine();
                        }
                        if(diepte==7){
                            if(args.length > 5){
                                if(lcaHeader.contains("root")){
                                    String[] root = new String[1];
                                    root[0] = "root";
                                    peptides.add(new Peptide(lcaHeader, taxonomy_score, root));
                                }else{
                                    String[] lineage = lineageR.readNext();
                                    peptides.add(new Peptide(lcaHeader, taxonomy_score, lineage));
                                }
                            }else{
                                peptides.add(new Peptide(lcaHeader,taxonomy_score));
                            }
                            printProteins((double)24/(double)frame.length(),proteins);
                        }
                        if(args.length > 5){
                            drawScoredFrame(frame,diepte+6,peptides,(diepte<5),trueLineage);
                        }
                        drawScoredFrame(frame,diepte,peptides,(diepte<5));
                        drawSplitPoints(frame,diepte,(diepte<5));
                        diepte+=1;
                    }else{
    //      Specifiek voor K-meren
                        ArrayList<Kmer> kmers = new ArrayList<>();
                        while(lca.hasNextLine() && lcaHeader.equals(header)){
                            String nextLCA = lca.nextLine();
                            if(nextLCA.startsWith(">")){
                                lcaHeader = nextLCA;
                            }else{
                                if(args.length > 5){
                                    if(nextLCA.contains("root")){
                                        String[] root = new String[1];
                                        root[0] = "root";
                                        kmers.add(new Kmer(nextLCA,k,taxonomy_score, root));;
                                    }else{
                                        String[] lineage = lineageR.readNext();
                                        kmers.add(new Kmer(nextLCA,k,taxonomy_score, lineage));;
                                    }
                                }else{
                                    kmers.add(new Kmer(nextLCA,k,taxonomy_score));
                                }
                            }
                        }
                        int[] frameDepths = getCoverageDepth(frame,kmers,k);
                        if(args.length > 5){
                            plotCoverageDepths(forward,diepte,frameDepths,k,2);
                            int framenr;
                            if((forward < 4)){  
                                framenr = forward;
                            }else{
                                framenr = 3-forward;
                            }
                            System.out.println("\\node[draw=none] at (24.5,-" + (double) (diepte/(double)2) + ") {" + framenr + "};");
                            diepte +=1;
                            drawKmerFrame(frame, diepte , (forward < 4),kmers,k);
                            diepte++;
                            drawKmerFrame(frame, diepte, (forward < 4),kmers, k, trueLineage);
                            diepte++;
                        }else{
                            plotCoverageDepths(forward,diepte,frameDepths,k,4);
                            diepte +=1;
                            diepte += drawKmerFrame(frame, diepte, forward ,kmers,k);
                        }
                        if(forward == 1){
                            printProteins((double)24/(double)frame.length(),proteins);
                        }
                        forward ++;
                    }
                }
            } catch (FileNotFoundException ex) {
                Logger.getLogger(ScoreReads.class.getName()).log(Level.SEVERE, null, ex);
            } catch (IOException ex) {
                Logger.getLogger(ScoreReads.class.getName()).log(Level.SEVERE, null, ex);
            }
            printFooter();
        }
    }
    
    public static void drawScoredFrame(String frame, int framenr, List<Peptide> peptides, boolean forward){
        double unit = (double) 24/(double)frame.length();
        for(Peptide pept:peptides){
            double start;
            double end;
            if(forward){
                start = (double) frame.indexOf(pept.aminoSeq)*unit;
                end = start + (double) pept.length*unit;
            }else{
                end = (double)frame.length()*unit-(double)frame.indexOf(pept.aminoSeq)*unit;
                start = end - (double) pept.length*unit;    
            }
            String taxonName = pept.taxonName;
            if(taxonName.equalsIgnoreCase("root")){
                System.out.println("\\draw[root] ("+start+",-"+framenr+") -- ("+end+",-"+framenr+");");
                double score = pept.getScore();
                System.out.println("\\node [below,font=\\tiny] at ("+(start+end)/(double) 2+",-"+framenr+") {"+score+"};");
            }else{
                System.out.println("\\draw["+pept.taxonRank+"] ("+start+",-"+framenr+") -- ("+end+",-"+framenr+");");
                double score = pept.getScore();
                System.out.println("\\node [below,font=\\tiny] at ("+(start+end)/(double) 2+",-"+framenr+") {"+score+"};");
            }
        }
    }
    
    public static void drawScoredFrame(String frame, int framenr, List<Peptide> peptides, boolean forward, String[] trueLineage){
        double unit = (double) 24/(double)frame.length();
        String lineageString="";
        for(int i=0; i<trueLineage.length;i++){
            lineageString=lineageString + trueLineage[i];
        }
        for (Peptide pept:peptides){
            double start;
            double end;
            if(forward){
                start = (double) frame.indexOf(pept.aminoSeq)*unit;
                end = start + (double) pept.length*unit;
            }else{
                end = (double)frame.length()*unit-(double)frame.indexOf(pept.aminoSeq)*unit;
                start = end - (double) pept.length*unit;    
            }
            String taxonName = pept.taxonName;
            if(taxonName.equalsIgnoreCase("root")){
                System.out.println("\\draw[root] ("+start+",-"+framenr+") -- ("+end+",-"+framenr+");");
            }else{
                if(lineageString.contains(taxonName)){
                    System.out.println("\\draw["+pept.taxonRank+"] ("+start+",-"+framenr+") -- ("+end+",-"+framenr+");");
                }else{
                    String[] lineage = pept.lineage;
                    printDisagreement(lineage,trueLineage,start,end,framenr);
                }
            }
        }
    }
    
    public static void drawKmerFrame(String frame, int framenr, boolean forward, List<Kmer> kmers, int k){
        double unit = (double) 24/(double)frame.length();
        int pos=0;
        for(Kmer kmer: kmers){
            double start;
            start = frame.indexOf(kmer.aminoSeq, pos);
            pos = (int) start + 1;
            start = (double) start*unit;
            if(! forward){
                start = (double)frame.length()*unit - start;   
            }
            String taxonName = kmer.taxonName;
            double above = ((double)framenr + 0.2)/2;
            double below = ((double)framenr - 0.2)/2;
            if(taxonName.equalsIgnoreCase("root")){
                System.out.println("\\draw[root, line cap=butt, line width = 2pt] ("+start+",-"+above+") -- ("+start+",-"+below+");");
            }else{
                System.out.println("\\draw["+kmer.taxonRank+", line cap=butt, line width = 2pt] ("+start+",-"+above+") -- ("+start+",-"+below+");");
            }
        }
    }
    public static void drawKmerFrame(String frame, int framenr, boolean forward, List<Kmer> kmers, int k, String[] trueLineage){
        String lineageString="";
        for(int i=0; i<trueLineage.length;i++){
            lineageString=lineageString + trueLineage[i];
        }
        double unit = (double) 24/(double)frame.length();
        int pos=0;
        for(Kmer kmer: kmers){
            double start;
            start = frame.indexOf(kmer.aminoSeq, pos);
            pos = (int) start + 1;
            start = (double) start*unit;
            if(! forward){
                start = (double)frame.length()*unit - start;   
            }
            String taxonName = kmer.taxonName;
            double above = ((double)framenr + 0.2)/2;
            double below = ((double)framenr - 0.2)/2;
            if(taxonName.equalsIgnoreCase("root")){
                System.out.println("\\draw[root, line cap=butt, line width = 2pt] ("+start+",-"+above+") -- ("+start+",-"+below+");");
            }else{
                if(lineageString.contains(taxonName)){
                    System.out.println("\\draw["+kmer.taxonRank+", line cap=butt, line width = 2pt] ("+start+",-"+above+") -- ("+start+",-"+below+");");
                }else{
                    String[] lineage = kmer.lineage;
                    printKmerDisagreement(lineage,trueLineage,start,framenr);
                }
            }
        }
    }
    
    public static int drawKmerFrame(String frame, int diepte, int forward, List<Kmer> kmers, int k){
        int depth = 0;
        double unit = (double) 24/(double)frame.length();
        int pos=0;
        int[] prevEnd = new int[k+1];
        int tekendiepte;
        for(Kmer kmer: kmers){
            tekendiepte = 0;
            double start;
            double end;
            start = frame.indexOf(kmer.aminoSeq, pos);
            pos = (int) start+1;
            while(tekendiepte < k && start<=prevEnd[tekendiepte]){
                tekendiepte ++;
            }
            prevEnd[tekendiepte] = (int) start + k -1;
            start = (double) start*unit;
            end = start + (double) (k-1)*unit;
            if(! (forward < 4)){
                end = (double)frame.length()*unit - start;
                start = end - (double) (k-1)*unit;    
            }
            String taxonName = kmer.taxonName;
            double framenr = (double) (diepte + tekendiepte)/4;
            if(taxonName.equalsIgnoreCase("root")){
                System.out.println("\\draw[root] ("+start+",-"+framenr+") -- ("+end+",-"+framenr+");");
            }else{
                System.out.println("\\draw["+kmer.taxonRank+"] ("+start+",-"+framenr+") -- ("+end+",-"+framenr+");");
            }
        }
        while(depth < k && prevEnd[depth]!=0){
            depth ++;
        }
        
        System.out.println("\\draw[frame] (24.2,-" + (double) (diepte-1.5)/4 + ") -- (24.2,-" + (double) (diepte + depth + 0.5)/4 + ");");
        System.out.println("\\draw[frame] (24,-" + (double) (diepte-1.5)/4 + ") -- (24.2,-" + (double) (diepte-1.5)/4+ ");");
        if(forward == 6){
            System.out.println("\\draw[frame] (24,-" + (double) (diepte + depth + 0.5)/4 + ") -- (24.2,-" + (double) (diepte + depth+0.5)/4 + ");");
        }
        int framenr;
        if((forward < 4)){  
            framenr = forward;
        }else{
            framenr = 3-forward;
        }
        System.out.println("\\node[draw=none] at (24.5,-" + (double) (diepte-1)/4 + ") {" + framenr + "};");
        return depth;
    }
    
    public static void drawSplitPoints(String frame, int framenr, boolean forward){
        double unit = (double) 24/(double)frame.length();
        double above = framenr + 0.2;
        double below = framenr - 0.2;
        int n = frame.length();
        int previous = 0;
        if(forward){
            for(int i = 1; i<n;i++){
               char ch = frame.charAt(i-1);
               boolean tryptic = (ch=='R'||ch=='K') && frame.charAt(i)!='P';
               if(ch == '*'||tryptic){
                    int length = i-previous;
                    if(length>4){
                        System.out.println("\\node [above,font=\\scriptsize] at ("+(double)(i+previous)*unit/(double) 2 + ",-"+framenr+") {"+length+"};");
                    }
                    previous = i;
                    System.out.println("\\draw [red] ("+i*unit+",-"+above+") -- ("+i*unit+",-"+below+");");
               }
            }
            int length = n-previous;
            if(length>4){
                System.out.println("\\node [above,font=\\scriptsize] at ("+(double)(n+previous)*unit/(double) 2 + ",-"+framenr+") {"+length+"};");
            }
        }else{
            for(int i = 1;i<n;i++){
                char ch = frame.charAt(i-1);
                boolean tryptic = (ch=='R'||ch=='K') && frame.charAt(i)!='P';
                if(ch == '*'||tryptic){
                    int length = i-previous;
                    if(length>4){
                        System.out.println("\\node [above,font=\\scriptsize] at ("+(n*unit - (double)(i+previous)*unit/(double) 2 )+ ",-"+framenr+") {"+length+"};");
                    }
                    previous = i;
                    System.out.println("\\draw [red] ("+(n-i)*unit+",-"+above+") -- ("+(n-i)*unit+",-"+below+");");
                }
            }
            int length = n-previous;
            if(length>4){
                System.out.println("\\node [above,font=\\scriptsize] at ("+(n*unit - (double)(n+previous)*unit/(double) 2 )+ ",-"+framenr+") {"+length+"};");
            }
        }
    }
    
    private static void fillTaxScore(String file) throws FileNotFoundException{
        taxonomy_score.put("no rank", 0.1);
        File input = new File(file);
        Scanner sc = new Scanner(input);
        while(sc.hasNextLine()){
            String line = "";
            String[] taxa = new String[0];
            while(!sc.hasNextDouble()){
                line = line + sc.next() + " ";
            }
            taxa = line.trim().split(",");
            double s = sc.nextDouble();
            for(String taxon : taxa){
                taxonomy_score.put(taxon, s);
            }
        }
    }
    
    
    private static int[] getCoverageDepth(String frame, List<Kmer> kmers, int k){
        int[] depths = new int[frame.length()];
        int cur_pos = 0;
        for(Kmer kmer:kmers){
            int start = frame.indexOf(kmer.aminoSeq, cur_pos);
            cur_pos = start;
            for(int i=start;i<start+k;i++){
                depths[i] = depths[i] + 1;
            }
        }
        return depths;
    }
    
    public static void plotCoverageDepths(int frame, int diepte, int[] depths, int k, int factor){
        double unit = (double) 24 / (double) depths.length;
        int n = depths.length;
        for(int i=0;i<depths.length;i++){
            int fact = depths[i]*100/k;
            double above = ((double)diepte + 0.2)/factor;
            double below = ((double)diepte - 0.2)/factor;
            if(frame<4){
                System.out.println("\\draw [RoyalBlue!"+fact+",line width=2pt] ("+(i)*unit+",-"+above+") -- ("+(i)*unit+",-"+below+");");
            }else{
                System.out.println("\\draw [RoyalBlue!"+fact+",line width=2pt] ("+(n-i)*unit+",-"+above+") -- ("+(n-i)*unit+",-"+below+");");
            }
        }
    }
    
    public static void printHeader(){
        System.out.println("\\begin{tikzpicture}\n" +
            "[font=\\sffamily, \n"+ 
            "protein/.style={violet, line width = 6pt, line cap = round},\n" +
            "frame/.style={orange, line width = 2pt, line cap = round},\n" +
            "root/.style={teal!10, line width = 4pt, line cap = round},\n" +
            "superkingdom/.style={teal!15, line width = 4pt, line cap = round},\n" +
            "kingdom/.style={teal!20, line width = 4pt, line cap = round},\n" +
            "subkingdom/.style={teal!25, line width = 4pt, line cap = round},\n" +
            "superphylum/.style={teal!30, line width = 4pt, line cap = round},\n" +
            "phylum/.style={teal!35, line width = 4pt, line cap = round},\n" +
            "subphylum/.style={teal!40, line width = 4pt, line cap = round},\n" +
            "superclass/.style={teal!45, line width = 4pt, line cap = round},\n" +
            "class/.style={teal!50, line width = 4pt, line cap = round},\n" +
            "subclass/.style={teal!55, line width = 4pt, line cap = round},\n" +
            "infraclass/.style={teal!60, line width = 4pt, line cap = round},\n" +
            "superorder/.style={teal!65, line width = 4pt, line cap = round},\n" +
            "order/.style={teal!70, line width = 4pt, line cap = round},\n" +
            "suborder/.style={teal!75, line width = 4pt, line cap = round},\n" +
            "infraorder/.style={teal!80, line width = 4pt, line cap = round},\n" +
            "parvorder/.style={teal!85, line width = 4pt, line cap = round},\n" +
            "superfamily/.style={teal!90, line width = 4pt, line cap = round},\n" +
            "family/.style={teal!95, line width = 4pt, line cap = round},\n" +
            "subfamily/.style={teal, line width = 4pt, line cap = round},\n" +
            "tribe/.style={teal!95!black, line width = 4pt, line cap = round},\n" +
            "subtribe/.style={teal!90!black, line width = 4pt, line cap = round},\n" +
            "genus/.style={teal!85!black, line width = 4pt, line cap = round},\n" +
            "subgenus/.style={teal!80!black, line width = 4pt, line cap = round},\n" +
            "species group/.style={teal!75!black, line width = 4pt, line cap = round},\n" +
            "species subgroup/.style={teal!70!black, line width = 4pt, line cap = round},\n" +
            "species/.style={teal!65!black, line width = 4pt, line cap = round},\n" +
            "subspecies/.style={teal!60!black, line width = 4pt, line cap = round},\n" +
            "varietas/.style={teal!55!black, line width = 4pt, line cap = round},\n" +
            "forma/.style={teal!50!black, line width = 4pt, line cap = round}]\n" +
            "\\node[font=\\bfseries\\LARGE,align=center,above] at (12,2) {"+title+"} ;\n" +
            "\\node[font=\\bfseries,align=center,above] at (12,1) {"+subtitle+"} ;\n" +
            "\\draw (0,0) -- (24,0) ;");
    }
    
    public static void printEmptyLines(int startindex){
        System.out.println(
        "\\draw (0,-"+startindex+") -- (24,-"+startindex+") node[anchor=south] {$+1$};\n" +
        "\\draw (0,-"+(startindex+1)+") -- (24,-"+(startindex+1)+") node[anchor=south] {$+2$};\n" +
        "\\draw (0,-"+(startindex+2)+") -- (24,-"+(startindex+2)+") node[anchor=south] {$+3$};\n" +
        "\\draw (0,-"+(startindex+3)+") -- (24,-"+(startindex+3)+") node[anchor=south] {$-1$};\n" +
        "\\draw (0,-"+(startindex+4)+") -- (24,-"+(startindex+4)+") node[anchor=south] {$-2$};\n" +
        "\\draw (0,-"+(startindex+5)+") -- (24,-"+(startindex+5)+") node[anchor=south] {$-3$};");
    }
    
    public static void printFooter(){
        System.out.println("\\end{tikzpicture}");
    }
    
    public static void printDisagreement(String[] lineage, String[] trueLin, double start, double end, int framenr){
        int[] disagreement = new int[2];
        int agreement = 0;
        while(agreement< lineage.length && agreement< trueLin.length &&lineage[agreement].trim().equals(trueLin[agreement])){
            agreement +=1;
        }
        disagreement[0] = lineage.length + 1 + trueLin.length - 2*agreement;
        disagreement[1] = 100-disagreement[0];
        System.out.println("\\draw[red!"+disagreement[1]+", line width = 4pt, line cap = round] ("+start+",-"+framenr+") -- ("+end+",-"+framenr+");");
        System.out.println("\\node[below, font=\\small] at ("+(start+end)/(double)2+",-"+framenr+") {"+disagreement[0]+"};");
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
        System.out.println("\\draw[red!"+disagreement[1]+", line width = 2pt, line cap = butt] ("+start+",-"+above+") -- ("+start+",-"+below+");");
        System.out.println("\\node[draw=none, below, font=\\sffamily\\fontsize{2}{1}\\selectfont] at ("+start+",-"+(below+0.2)+") {"+disagreement[0]+"};");
    }
    
    private static void printProteins(double unit, String[] proteinPos) {
        if(proteinPos.length!=0){
            for (int i = 0; i < proteinPos.length-1;i+=2){
                double proteinstart = Double.parseDouble(proteinPos[i])*unit / (double) 3;
                double proteinstop = Double.parseDouble(proteinPos[i+1])*unit / (double) 3;
                System.out.println("\\draw[protein] ("+proteinstart+",0) -- ("+proteinstop+",0);");
            }
        }
    }
}
