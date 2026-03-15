package at.univie.centralized;

import at.univie.enumhyp.GenerateEnumerate;

import java.io.IOException;

/**
 * Centralized approach - processes data without distribution.
 * Uses centralized/data directory for input/output.
 */
public class CentralizedBaseline {
    
//    public static void main(String[] args) throws IOException, InterruptedException {
//
//
//        // Create GenerateEnumerate instance
//        GenerateEnumerate genEnum = new GenerateEnumerate();
//
//        // Process flight data
//        String csvFile = "fdReduced";
//        boolean genGraph = true;
//        boolean genTxt = true;
//
//        genEnum.generateEnumerateSingle(csvFile, genGraph, genTxt);
//
//        System.out.println("Centralized processing completed for: " + csvFile);
//        System.out.println("Results are in: centralized/data/" + csvFile + "/");
//    }

    public static void main(String[] args) throws IOException, InterruptedException {

        // Create GenerateEnumerate instance
        GenerateEnumerate genEnum = new GenerateEnumerate();

        // Process flight data
        String csvFile = "abalone";
        boolean genGraph = true;
        boolean genTxt = true;

        // ---- TIMER START ----
        long startNs = System.nanoTime();

        genEnum.generateEnumerateSingle(csvFile, genGraph, genTxt);

        // ---- TIMER END ----
        long endNs = System.nanoTime();

        long elapsedMs = (endNs - startNs) / 1_000_000L;
        double elapsedSec = (endNs - startNs) / 1_000_000_000.0;

        System.out.println("[TIMER] generateEnumerateSingle(" + csvFile + "): " + elapsedMs + " ms (" + elapsedSec + " s)");

        System.out.println("Centralized processing completed for: " + csvFile);
        System.out.println("Results are in: centralized/data/" + csvFile + "/");
    }

}
