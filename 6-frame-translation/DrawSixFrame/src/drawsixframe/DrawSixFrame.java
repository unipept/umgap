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
    
    
    public static void organismUnknownPrint(String lcaFile, String sixframeFile){
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
        "\\draw node[anchor=north] {True sequence} (0,0) -- (24,0) ;\n" +
        "\\draw (0,-2) -- (24,-2) node[anchor=south] {$+1$};\n" +
        "\\draw (0,-3) -- (24,-3) node[anchor=south] {$+2$};\n" +
        "\\draw (0,-4) -- (24,-4) node[anchor=south] {$+3$};\n" +
        "\\draw (0,-5) -- (24,-5) node[anchor=south] {$-1$};\n" +
        "\\draw (0,-6) -- (24,-6) node[anchor=south] {$-2$};\n" +
        "\\draw (0,-7) -- (24,-7) node[anchor=south] {$-3$};");
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
            System.out.println("\\end{tikzpicture}");
            } catch (IOException ex) {
            System.out.println("File not found");
        }
    }
    
    public static void organismKnownPrint(String lcaFile, String sixFrameFile, String lineageFile){
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
        "\\draw node[anchor=north] {True sequence} (0,0) -- (24,0) ;\n" +
        "\\draw (0,-2) -- (24,-2) node[anchor=south] {$+1$};\n" +
        "\\draw (0,-3) -- (24,-3) node[anchor=south] {$+2$};\n" +
        "\\draw (0,-4) -- (24,-4) node[anchor=south] {$+3$};\n" +
        "\\draw (0,-5) -- (24,-5) node[anchor=south] {$-1$};\n" +
        "\\draw (0,-6) -- (24,-6) node[anchor=south] {$-2$};\n" +
        "\\draw (0,-7) -- (24,-7) node[anchor=south] {$-3$};");
        String[] trueLin = correctTaxonomy.split("; ");
        CSVReader reader1;
        CSVReader reader2;
        try{
            reader1 = new CSVReader(new FileReader(lineageFile),';');
            List<String[]> lineages = reader1.readAll();
            Iterator<String[]> linIt = lineages.iterator();
            reader2 = new CSVReader(new FileReader(lcaFile));
            BufferedReader sixframe = new BufferedReader(new FileReader(sixFrameFile));
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
                                linIt.next();
                                System.out.println("\\draw["+nextline[4]+"] ("+start+",-"+framenr+") -- ("+end+",-"+framenr+");");
                            }else{
                                String[] lineage = linIt.next();
                                int disagreement = getDisagreement(lineage,trueLin);
                                System.out.println("\\draw[red!"+disagreement+", line width = 4pt, line cap = round] ("+start+",-"+framenr+") -- ("+end+",-"+framenr+");");
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
                                int disagreement = getDisagreement(lineage,trueLin);
                                System.out.println("\\draw[red!"+disagreement+", line width = 4pt, line cap = round] ("+start+",-"+framenr+") -- ("+end+",-"+framenr+");");
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
                    int disagreement = getDisagreement(lineage,trueLin);
                    System.out.println("\\draw[red!"+disagreement+", line width = 4pt, line cap = round] ("+start+",-"+framenr+") -- ("+end+",-"+framenr+");");
                }
            }
            System.out.println("\\end{tikzpicture}");
        }catch(Exception ex){
            ex.printStackTrace();
        }
    }
    
    public static int getDisagreement(String[] lineage, String[] trueLin){
        int agreement = 0;
        while(agreement< lineage.length && agreement< trueLin.length &&lineage[agreement].trim().equals(trueLin[agreement])){
            agreement +=1;
        }
        int disagreement = lineage.length + trueLin.length - 2*agreement;
        return 100 - disagreement*2;
    }
    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
//        System.out.println("\\begin{tikzpicture}\n" +
//        "[protein/.style={violet, line width = 6pt, line cap = round},\n" +
//        "root/.style={cyan, line width = 4pt, line cap = round},\n" +
//        "genus/.style={teal!60, line width = 4pt, line cap = round},\n" +
//        "species/.style={teal, line width = 4pt, line cap = round},\n" +
//        "species group/.style={teal!80, line width = 4pt, line cap = round},\n" +
//        "wrong/.style={pink, line width = 4pt, line cap = round},\n" +
//        "superkingdom/.style={teal!10, line width = 4pt, line cap = round}]\n" +
//        "\\draw node[anchor=north] {True sequence} (0,0) -- (24,0) ;\n" +
//        "\\draw (0,-2) -- (24,-2) node[anchor=south] {$+1$};\n" +
//        "\\draw (0,-3) -- (24,-3) node[anchor=south] {$+2$};\n" +
//        "\\draw (0,-4) -- (24,-4) node[anchor=south] {$+3$};\n" +
//        "\\draw (0,-5) -- (24,-5) node[anchor=south] {$-1$};\n" +
//        "\\draw (0,-6) -- (24,-6) node[anchor=south] {$-2$};\n" +
//        "\\draw (0,-7) -- (24,-7) node[anchor=south] {$-3$};");
        proteinSeq.add("MKISDLMTYHGCKNRKELSEKTGYSTVTLWKWENNGIPARTQAVLQVKTKGKLKADLQALTA");
        proteinSeq.add("MSLHSRIRQKLEEKKLRAADLARATKKSPVAVKKWLDGTSVPTAENLKVIAKFLGVSDDWLLYGGPVEQESNNLPQLNVLDIEAFKQKYNIPDSEDAVKFVQTSDKPFPIQKRYVPVKAYSKMGMDGYFTDMGYEGNAGDGYVPTHSAGPRAYGIKGTGDSMFPAIRNGWYVVCDPDAELVPTEFVQVCLKDGRCTIKEFIGINNDVLSLIAVNGGERLSFNMDEVESITAITDIVPPSQHRQEHPYSH");
        taxonomy.addAll(Arrays.asList(correctTaxonomy.split("; ")));
//        organismUnknownPrint("found_lcas_double.txt","double.sixframe");
        organismKnownPrint("found_lcas_double.txt","double.sixframe","double.lineage");
//        CSVReader reader;
//        try {
//            reader = new CSVReader(new FileReader("found_lcas_double.txt"));
//            BufferedReader sixframe = new BufferedReader(new FileReader("double.sixframe"));
//            List<String[]> lcas = reader.readAll();
//            Iterator<String[]> it = lcas.iterator();
//            it.next();
//            String header;
//            int framenr = 1;
//            String[] nextline = it.next();
//            String frame = "";
//            double unit = 1;
//            while((header=sixframe.readLine())!=null){
//                framenr+=1;
//                frame = sixframe.readLine();
//                unit = pagewidth/(double)frame.length();
//                for(String seq: proteinSeq){
//                    if(frame.contains(seq)){
//                        double proteinstart;
//                        double proteinstop;
//                        if(framenr>=2 && framenr <5){
//                            proteinstart = (double) frame.indexOf(seq)*unit;
//                            proteinstop = proteinstart + (double) seq.length()*unit;
//                        }else{
//                            proteinstop = (double)frame.length()*unit-(double)frame.indexOf(seq)*unit;
//                            proteinstart = proteinstop - (double) seq.length()*unit;
//                        }
//                    System.out.println("\\draw[protein] ("+proteinstart+",0) -- ("+proteinstop+",0);");
//                    }
//                }
//                boolean sameHeader = nextline[0].equals(header);
//                while(it.hasNext() && sameHeader){
//                    if(framenr>=2 && framenr <5){
//                        double start = (double) frame.indexOf(nextline[1])*unit;
//                        double end = start + (double) nextline[1].length()*unit;
//                        String taxonName = nextline[3];
//                        if(taxonName.equalsIgnoreCase("root")){
//                            System.out.println("\\draw[root] ("+start+",-"+framenr+") -- ("+end+",-"+framenr+");");
//                        }else{
//                            if(taxonomy.contains(taxonName)){
//                                System.out.println("\\draw["+nextline[4]+"] ("+start+",-"+framenr+") -- ("+end+",-"+framenr+");");
//                            }else{
//                                System.out.println("\\draw[wrong] ("+start+",-"+framenr+") -- ("+end+",-"+framenr+");");
//                            }
//                        }
//                    }else{
//                        double end = (double)frame.length()*unit-(double)frame.indexOf(nextline[1])*unit;
//                        double start = end - (double) nextline[1].length()*unit;
//                        String taxonName = nextline[3];
//                        if(taxonName.equalsIgnoreCase("root")){
//                            System.out.println("\\draw[root] ("+start+",-"+framenr+") -- ("+end+",-"+framenr+");");
//                        }else{
//                            if(taxonomy.contains(taxonName)){
//                                System.out.println("\\draw["+nextline[4]+"] ("+start+",-"+framenr+") -- ("+end+",-"+framenr+");");
//                            }else{
//                                System.out.println("\\draw[wrong] ("+start+",-"+framenr+") -- ("+end+",-"+framenr+");");
//                            }
//                        }
//                    }
//                    nextline = it.next();
//                    sameHeader = nextline[0].equals(header);
//                }
//            }
//            double end = (double)frame.length()*unit-(double)frame.indexOf(nextline[1])*unit;
//            double start = end - (double) nextline[1].length()*unit;
//            String taxonName = nextline[3];
//            if(taxonName.equalsIgnoreCase("root")){
//                System.out.println("\\draw[root] ("+start+",-"+framenr+") -- ("+end+",-"+framenr+");");
//            }else{
//                if(taxonomy.contains(taxonName)){
//                    System.out.println("\\draw["+nextline[4]+"] ("+start+",-"+framenr+") -- ("+end+",-"+framenr+");");
//                }else{
//                    System.out.println("\\draw[wrong] ("+start+",-"+framenr+") -- ("+end+",-"+framenr+");");
//                }
//            }
//            System.out.println("\\end{tikzpicture}");
//            } catch (IOException ex) {
//            System.out.println("File found_lcas.txt not found");
//        }
    }
    
}
