/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package seedextend;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Deque;
import java.util.Iterator;
import java.util.List;
import java.util.Scanner;
import java.util.TreeMap;

/**
 *
 * @author Aranka
 */
public class SeedExtend {
    private static final int minSeedSize = 3;
    private static int k;
    private static List<ExtendedSeed> extendedSeedsPerFrame;

    /**
     * @param args the command line arguments
     * 0: =0 voor peptiden =k voor k-meren
     * 1: sixframe-file
     * 2: lca-file
     *              
     */
    public static void main(String[] args) throws IOException {
        BufferedReader sixframe = null;
        try {
            k = Integer.parseInt(args[0]);
            File sixframeF = new File(args[1]);
            sixframe = new BufferedReader(new FileReader(sixframeF));
            Scanner lcaF = new Scanner(new File(args[2]));
            String lcaHeader = lcaF.nextLine();
            while(! lcaHeader.isEmpty()){
                String printHeader = lcaHeader.substring(0,lcaHeader.indexOf("|"));
                int frameN = 1;
                extendedSeedsPerFrame = new ArrayList<>();
                while(frameN <= 6){
                    Frame f = new Frame(frameN,k);
                    String header;
                    // Initialiseer lijst met seeds, index en deque van kmeren en lees de eerste regels van beide files
                    TreeMap<Integer,Seed> frameSeeds = new TreeMap<>();
                    TreeMap<Integer,Kmer> frameKmers = new TreeMap<>();
                    int prevIndex=0;
                    header = sixframe.readLine();
                    String frame = sixframe.readLine();
                    Deque<Kmer> deq = new ArrayDeque<>();
                    while(lcaHeader.equals(header)&&lcaF.hasNextLine()){
                        // voor alle lca's bij deze frame
                        String nextLCA = lcaF.nextLine();
                        if(nextLCA.startsWith(">")){
                            lcaHeader = nextLCA;
                            if(deq.size() >= 3){
                                Seed newSeed = new Seed(deq);
                                frameSeeds.put(newSeed.start,newSeed);
                            }
                        }else{
                            // maak de nieuwe kmer
                            Kmer kmer = new Kmer(nextLCA,k);
                            int start = frame.indexOf(kmer.aminoSeq, prevIndex);
                            kmer.setStart(start);
                            frameKmers.put(start,kmer);
                            prevIndex = start;
                            if (!deq.isEmpty()){
                                // voorwaarde om bij de huidige seed te mogen horen
                                if (deq.peekLast().start == start - 1 && deq.peekLast().taxonID == kmer.taxonID){
                                    deq.add(kmer);
                                }else{
                                    // minimale grootte voor een seed
                                    if(deq.size() >= 3){
                                        Seed newSeed = new Seed(deq);
                                        frameSeeds.put(newSeed.start,newSeed);
                                    }
                                    deq.clear();
                                    deq.add(kmer);
                                }
                            }else{
                                deq.add(kmer);
                            }
                        }
                    }
                    f.fillKmers(frameKmers);
                    f.fillSeeds(frameSeeds);
                    f.extendSeeds();
                    ArrayList<ExtendedSeed> frameResult = f.getExtendedSeeds();
                    extendedSeedsPerFrame.addAll(frameResult);
                    frameN ++;
                }
                if(! extendedSeedsPerFrame.isEmpty()){
                    printBestSeed(printHeader);
                }
            }
        } catch (NumberFormatException | IOException ex) {
            System.out.println("One of both files could not be read");
        }
        //System.out.println("Done");
    }
    
    private static void printBestSeed(String header){
        ExtendedSeed longest = extendedSeedsPerFrame.get(0);
        for(int i=1; i< extendedSeedsPerFrame.size();i++){
            if(extendedSeedsPerFrame.get(i).getLength()>longest.getLength()){
                longest = extendedSeedsPerFrame.get(i);
            }
        }
        System.out.println(header + "|Frame " + longest.frameN);
        System.out.println(longest.start + "," + longest.end + ": " + longest.taxonID);
    }
}
