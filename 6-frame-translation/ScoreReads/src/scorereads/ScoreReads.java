/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package scorereads;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
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
     */
    public static void main(String[] args) {
        printHeader();
        try {
            fillTaxScore("taxonomy_score.txt");
        } catch (FileNotFoundException ex) {
            Logger.getLogger(ScoreReads.class.getName()).log(Level.SEVERE, null, ex);
        }
        File sixFrame = new File("single_protein.sixframe");
        try {
            BufferedReader sixframe = new BufferedReader(new FileReader(sixFrame));
            Scanner lca = new Scanner(new File("lca_file.temp"));
            String lcaHeader = lca.nextLine();
//            String unipeptInfo = lca.nextLine();
//            int framecount = 2;
            int diepte = 2;
            int forward = 1;
            String header;
            int k=9;
            while((header=sixframe.readLine())!=null){
                String frame = sixframe.readLine();
//    Specifiek voor Triptische Peptiden            
//                ArrayList<Peptide> peptides = new ArrayList<>();
//                while(lca.hasNextLine() && unipeptInfo.contains(header)){
//                    peptides.add(new Peptide(unipeptInfo,taxonomy_score));
//                    unipeptInfo = lca.nextLine();
//                }
//                if(framecount==13){
//                    peptides.add(new Peptide(unipeptInfo,taxonomy_score));
//                }
//                drawScoredFrame(frame,framecount,peptides,(framecount<11));
//                drawSplitPoints(frame,framecount,(framecount<11));
//                framecount+=1;
//      Specifiek voor K-meren
                ArrayList<Kmer> kmers = new ArrayList<>();
                while(lca.hasNextLine() && lcaHeader.equals(header)){
                    String nextLCA = lca.nextLine();
                    if(nextLCA.startsWith(">")){
                        lcaHeader = nextLCA;
                    }else{
                        kmers.add(new Kmer(nextLCA,k,taxonomy_score));
                    }
                }
                diepte += drawKmerFrame(frame, diepte, forward ,kmers,k);
                forward ++;
            }
        } catch (FileNotFoundException ex) {
            Logger.getLogger(ScoreReads.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(ScoreReads.class.getName()).log(Level.SEVERE, null, ex);
        }
        printFooter();
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
        
        System.out.println("\\draw[frame] (24.2,-" + (double) (diepte-0.5)/4 + ") -- (24.2,-" + (double) (diepte + depth + 0.5)/4 + ");");
        System.out.println("\\draw[frame] (24,-" + (double) (diepte-0.5)/4 + ") -- (24.2,-" + (double) (diepte-0.5)/4+ ");");
        if(forward == 6){
            System.out.println("\\draw[frame] (24,-" + (double) (diepte + depth + 0.5)/4 + ") -- (24.2,-" + (double) (diepte + depth+0.5)/4 + ");");
        }
        int framenr;
        if((forward < 4)){  
            framenr = forward;
        }else{
            framenr = 3-forward;
        }
        System.out.println("\\node[draw=none] at (24.5,-" + (double) diepte/4 + ") {" + framenr + "};");
        return depth;
    }
    
    public static void drawSplitPoints(String frame, int framenr, boolean forward){
        double unit = (double) 24/(double)frame.length();
        double above = framenr + 0.2;
        double below = framenr - 0.2;
        int n = frame.length();
        int previous = 0;
        if(forward){
            for(int i = 1; i<=n;i++){
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
            for(int i = 1;i<=n;i++){
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
    public static void printHeader(){
        System.out.println("\\begin{tikzpicture}\n" +
        "[protein/.style={violet, line width = 6pt, line cap = round},\n" +
        "root/.style={teal!20, line width = 4pt, line cap = round},\n" +
        "genus/.style={teal, line width = 4pt, line cap = round},\n" +
        "varietas/.style={teal!50!black, line width = 4pt, line cap = round},\n" +
        "subspecies/.style={teal!60!black, line width = 4pt, line cap = round},\n" +
        "species/.style={teal!70!black, line width = 4pt, line cap = round},\n" +
        "species group/.style={teal!80!black, line width = 4pt, line cap = round},\n" +
        "family/.style={teal!80, line width = 4pt, line cap = round},\n" +
        "suborder/.style={teal!75, line width = 4pt, line cap = round},\n" +
        "order/.style={teal!70, line width = 4pt, line cap = round},\n" +
        "superorder/.style={teal!65, line width = 4pt, line cap = round},\n" +
        "class/.style={teal!60, line width = 4pt, line cap = round},\n" +
        "subphylum/.style={teal!55, line width = 4pt, line cap = round},\n" +
        "phylum/.style={teal!50, line width = 4pt, line cap = round},\n" +
        "subkingdom/.style={teal!45, line width = 4pt, line cap = round},\n" +
        "kingdom/.style={teal!40, line width = 4pt, line cap = round},\n" +
        "superkingdom/.style={teal!30, line width = 4pt, line cap = round}]\n" +
        "\\node[font=\\bfseries\\LARGE,align=center,above] at (12,2) {"+title+"} ;\n" +
        "\\node[font=\\bfseries,align=center,above] at (12,1) {"+subtitle+"} ;\n" +
        "\\draw (0,0) -- (24,0) ;\n" +
        "\\draw (0,-2) -- (24,-2) node[anchor=south] {$+1$};\n" +
        "\\draw (0,-3) -- (24,-3) node[anchor=south] {$+2$};\n" +
        "\\draw (0,-4) -- (24,-4) node[anchor=south] {$+3$};\n" +
        "\\draw (0,-5) -- (24,-5) node[anchor=south] {$-1$};\n" +
        "\\draw (0,-6) -- (24,-6) node[anchor=south] {$-2$};\n" +
        "\\draw (0,-7) -- (24,-7) node[anchor=south] {$-3$};");
    }
    
    public static void printFooter(){
        System.out.println("\\end{tikzpicture}");
    }
}
