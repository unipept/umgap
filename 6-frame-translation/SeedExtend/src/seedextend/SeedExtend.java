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
import java.util.List;
import java.util.Scanner;
import java.util.TreeMap;

/**
 *
 * @author Aranka
 */
public class SeedExtend {
    private static int minSeedSize;
    private static int gapSize;
    private static double gapPenalty;
    private static int k;
    private static final List<ExtendedSeed> extendedSeeds = new ArrayList<>();
    private static final ArrayList<Frame> frames = new ArrayList<>();

    /**
     * @param args the command line arguments
     * 0: =0 voor peptiden =k voor k-meren
     * 1: sixframe-file
     * 2: lca-file
     * 3: minimal seed size
     * 4: gap size
     * 5: gap penalty
     */
    public static void main(String[] args) throws IOException {
        BufferedReader sixframe = null;
        minSeedSize = Integer.parseInt(args[3]);
        gapSize = Integer.parseInt(args[4]);
        gapPenalty = Double.parseDouble(args[5]);
        try {
            k = Integer.parseInt(args[0]);
            File sixframeF = new File(args[1]);
            sixframe = new BufferedReader(new FileReader(sixframeF));
            Scanner lcaF = new Scanner(new File(args[2]));
            String lcaHeader = lcaF.nextLine();
            while(lcaF.hasNextLine()){
                String printHeader = lcaHeader.substring(0,lcaHeader.lastIndexOf("|"));
                int frameN = 1;
                while(frameN <= 6){
                    Frame f = new Frame(frameN,k,gapSize,gapPenalty);
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
                            if(deq.size() >= minSeedSize){
                                Seed newSeed = new Seed(deq, frameN);
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
                                    if(deq.size() >= minSeedSize){
                                        Seed newSeed = new Seed(deq, frameN);
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
                    frames.add(f);
                    frameN ++;
                }
                getExtendedSeeds();
//                Zoek hoogst scorende extended seed en print de aanwezige taxa uit (VOOR SCORE2)
//                double maxScore = 0;
//                ExtendedSeed best = null;
//                for(ExtendedSeed e:extendedSeeds){
//                    double escore = e.getScore();
//                    if(escore > maxScore){
//                        maxScore = escore;
//                        best = e;
//                    }
//                }
//                if(best != null){
//                    System.out.println(printHeader);
//                    for (Kmer kmer : best.kmers) {
//                        System.out.print(kmer.taxonID + " ");
//                    }
//                }
//                System.out.println();
//                METHODE VOOR SCORE1, alles uitprinten
                if(! extendedSeeds.isEmpty()){
                    System.out.println(printHeader);
                    for(ExtendedSeed e:extendedSeeds){
                        for (Kmer kmer : e.kmers) {
                            System.out.print(kmer.taxonID + " ");
                        }
                    }
                    System.out.println();
                }
                extendedSeeds.clear();
                frames.clear();
            }
        } catch (NumberFormatException | IOException ex) {
            System.out.println("One of both files could not be read");
        }
        
    }
    
    private static void getExtendedSeeds(){  
        while(! frames.get(0).getSeeds().isEmpty()
                || ! frames.get(1).getSeeds().isEmpty()
                || ! frames.get(2).getSeeds().isEmpty()
                || ! frames.get(3).getSeeds().isEmpty()
                || ! frames.get(4).getSeeds().isEmpty()
                || ! frames.get(5).getSeeds().isEmpty()){
            double bestScore = 0;
            int[] bestPos = new int[2];
            for(int i = 0; i < 6; i ++){
                TreeMap<Integer,Seed> seeds = frames.get(i).getSeeds();
                for(Integer pos: seeds.keySet()){
                    Seed test = seeds.get(pos);
                    if(bestScore < test.getScore()){
                        bestScore = test.getScore();
                        bestPos[0] = i;
                        bestPos[1] = pos;
                    }
                }
            }
            extendedSeeds.add(frames.get(bestPos[0]).extendSeed(bestPos[1]));
        }
    }
    
}
