/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package drawsixframe;

import com.opencsv.CSVReader;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Scanner;
import java.util.TreeSet;
import java.util.logging.Level;
import java.util.logging.Logger;


/**
 *
 * @author Aranka
 */
public class DrawSixFrame {
    private static final String correctTaxonomy = "cellular organisms; Bacteria; Proteobacteria; Gammaproteobacteria; Pseudomonadales; Moraxellaceae; Acinetobacter; Acinetobacter calcoaceticus/baumannii complex; Acinetobacter baumannii";
    private static final TreeSet<String> taxonomy = new TreeSet<>();
    private static final double pagewidth = (double) 24;
    private static final HashMap<String,String> proteinSeq = new HashMap<>();
    private static final String title = "Acinetobacter Baumannii";
    private static String subtitle;
    private static final String folder = "./acinetobacter/";
    
    
    public static void organismUnknownPrint(String lcaFile, String sixframeFile){
        CSVReader reader;
        try {
            reader = new CSVReader(new FileReader(lcaFile));
            BufferedReader sixframe = new BufferedReader(new FileReader(sixframeFile));
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
                printProteins(frame,unit,framenr);
                boolean sameHeader = nextline[0].equals(header);
                while(it.hasNext() && sameHeader){
                    if(framenr>=2 && framenr <5){
                        double start = (double) frame.indexOf(nextline[1])*unit;
                        double end = start + (double) nextline[1].length()*unit;
                        String taxonName = nextline[3];
                        if(taxonName.equalsIgnoreCase("root")){
                            System.out.println("\\draw[root] ("+start+",-"+framenr+") -- ("+end+",-"+framenr+");");
                        }else{
                            System.out.println("\\draw["+nextline[4]+"] ("+start+",-"+framenr+") -- ("+end+",-"+framenr+");");
//                            System.out.println("\\node [below,font=\\footnotesize] at ("+end+",-"+framenr+") {"+taxonName+"};");
                        }
                    }else{
                        double end = (double)frame.length()*unit-(double)frame.indexOf(nextline[1])*unit;
                        double start = end - (double) nextline[1].length()*unit;
                        String taxonName = nextline[3];
                        if(taxonName.equalsIgnoreCase("root")){
                            System.out.println("\\draw[root] ("+start+",-"+framenr+") -- ("+end+",-"+framenr+");");
                        }else{
                            System.out.println("\\draw["+nextline[4]+"] ("+start+",-"+framenr+") -- ("+end+",-"+framenr+");");
//                            System.out.println("\\node [below,font=\\footnotesize] at ("+end+",-"+framenr+") {"+taxonName+"};");
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
                System.out.println("\\draw["+nextline[4]+"] ("+start+",-"+framenr+") -- ("+end+",-"+framenr+");");
//                System.out.println("\\node [below,font=\\footnotesize] at ("+end+",-"+framenr+") {"+taxonName+"};");
            }
            } catch (IOException ex) {
            System.out.println("File not found");
        }
    }
    
    public static void printHeader(){
        System.out.println("\\begin{tikzpicture}\n" +
        "[protein/.style={violet, line width = 6pt, line cap = round},\n" +
        "root/.style={teal!20, line width = 4pt, line cap = round},\n" +
        "genus/.style={teal, line width = 4pt, line cap = round},\n" +
        "subspecies/.style={teal!90!black, line width = 4pt, line cap = round},\n" +
        "species/.style={teal!70!black, line width = 4pt, line cap = round},\n" +
        "species group/.style={teal!80!black, line width = 4pt, line cap = round},\n" +
        "family/.style={teal!80, line width = 4pt, line cap = round},\n" +
        "order/.style={teal!70, line width = 4pt, line cap = round},\n" +
        "class/.style={teal!60, line width = 4pt, line cap = round},\n" +
        "phylum/.style={teal!50, line width = 4pt, line cap = round},\n" +
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
    
    public static void organismKnownPrint(String lcaFile, String sixFrameFile, String lineageFile){
        String[] trueLin = correctTaxonomy.split("; ");
        CSVReader reader1;
        CSVReader reader2;
        try{
            File lineageF = new File(lineageFile);
            reader1 = new CSVReader(new FileReader(lineageF),';');
            List<String[]> lineages = reader1.readAll();
            Iterator<String[]> linIt = lineages.iterator();
            File lcaF = new File(lcaFile);
            reader2 = new CSVReader(new FileReader(lcaF));
            File sixFrameF = new File(sixFrameFile);
            BufferedReader sixframe = new BufferedReader(new FileReader(sixFrameF));
            List<String[]> lcas = reader2.readAll();
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
                printProteins(frame,unit,framenr);
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
                                linIt.next();
                                System.out.println("\\draw["+nextline[4]+"] ("+start+",-"+framenr+") -- ("+end+",-"+framenr+");");
                            }else{
                                String[] lineage = linIt.next();
                                getDisagreement(lineage,trueLin,start,end,framenr);
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
                                linIt.next();
                                System.out.println("\\draw["+nextline[4]+"] ("+start+",-"+framenr+") -- ("+end+",-"+framenr+");");
                            }else{
                                String[] lineage = linIt.next();
                                getDisagreement(lineage,trueLin,start,end,framenr);
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
                    linIt.next();
                    System.out.println("\\draw["+nextline[4]+"] ("+start+",-"+framenr+") -- ("+end+",-"+framenr+");");
                }else{
                    String[] lineage = linIt.next();
                    getDisagreement(lineage,trueLin, start, end, framenr);
                }
            }
        }catch(Exception ex){
            ex.printStackTrace();
        }
    }
    
    public static void getDisagreement(String[] lineage, String[] trueLin, double start, double end, int framenr){
        int[] disagreement = new int[2];
        int agreement = 0;
        while(agreement< lineage.length && agreement< trueLin.length &&lineage[agreement].trim().equals(trueLin[agreement])){
            agreement +=1;
        }
        disagreement[0] = lineage.length + 1 + trueLin.length - 2*agreement;
        disagreement[1] = 100-disagreement[0]*2;
        System.out.println("\\draw[red!"+disagreement[1]+", line width = 4pt, line cap = round] ("+start+",-"+framenr+") -- ("+end+",-"+framenr+");");
        System.out.println("\\node[below] at ("+end+",-"+framenr+") {"+disagreement[0]+"};");
    }
    /**
     * @param args the command line arguments
     */
    public static void main(String[] args){
        boolean organismKnown = false;
        subtitle = "Translation table: " + args[0];
        String foundLCAS =  args[1];
        String sixframe =  args[2];
        String[] proteinFiles = new String[0];
        if(!args[3].equals("none")){
            proteinFiles= args[3].split(",");
        }
        String lineage = "";
        if(args.length > 4){
            organismKnown = true;
            lineage = args[4];
        }
        for(int i=0;i<proteinFiles.length;i++){
            try {
                Scanner sc = new Scanner(new File(proteinFiles[i]));
                String header = sc.nextLine();
                String description = header.split(" ")[0];
                String sequence = "";
                while(sc.hasNextLine()){
                    sequence = sequence + sc.nextLine();
                }
                proteinSeq.put(sequence,description);
            } catch (FileNotFoundException ex) {
                Logger.getLogger(DrawSixFrame.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        taxonomy.addAll(Arrays.asList(correctTaxonomy.split("; ")));
        printHeader();
        if(!organismKnown){
            organismUnknownPrint(foundLCAS,sixframe);
        }else{
            organismKnownPrint(foundLCAS,sixframe,lineage);
        }
        printFooter();
    }

    public static void printProteins(String frame, double unit, int framenr) {
        for(String seq: proteinSeq.keySet()){
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
                        System.out.println("\\node[below] at ("+proteinstop+",0)  {"+proteinSeq.get(seq)+"};");
                    }
                }
    }
    
}
