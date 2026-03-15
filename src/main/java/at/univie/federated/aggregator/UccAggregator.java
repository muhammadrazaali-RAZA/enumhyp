package at.univie.federated.aggregator;

import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;

/**
 * Aggregates local UCC outputs into global UCCs.
 * 
 * Algorithm: Tr(H1 ∪ ... ∪ Hn) = min { T1 ∪ ... ∪ Tn | Ti ∈ Tr(Hi) }
 * Implemented iteratively with minimization after each node.
 *
 */
public final class UccAggregator {
    
    private final FilePathManager pathManager;
    private final UccFileReader fileReader;
    private final UccMinimizer minimizer;
    
    public UccAggregator(Path inputDir, String filePrefix, String fileSuffix) {
        this.pathManager = new FilePathManager(inputDir, filePrefix, fileSuffix);
        this.fileReader = new UccFileReader();
        this.minimizer = new UccMinimizer();
    }
    
    /**
     * Builds paths for known number of parts.
     */
    public List<Path> buildPartPaths(int n) {
        return pathManager.buildPartPaths(n);
    }
    
    /**
     * Auto-discovers part file paths.
     */
    public List<Path> discoverPartPaths() throws IOException {
        return pathManager.discoverPartPaths();
    }
    
    /**
     * Reads a local UCC file.
     */
    public UccFileReader.LocalUccResult readLocalUccFile(Path file) throws IOException {
        return fileReader.readLocalUccFile(file);
    }
    
    /**
     * Aggregates global UCCs from local UCC lists.
     * Starts with {∅}, then for each node: candidates = {g ∪ t | g in global, t in local_i}, then minimalize.
     */
    public List<BitSet> aggregateGlobalUccs(List<List<BitSet>> localUccLists, int attributeCount) {
        // Start with empty set
        List<BitSet> global = new ArrayList<>(1);
        global.add(new BitSet(attributeCount));
        
        for (List<BitSet> local : localUccLists) {
            List<BitSet> candidates = buildCandidates(global, local);
            global = minimizer.minimalize(candidates);
        }
        
        return global;
    }
    
    private List<BitSet> buildCandidates(List<BitSet> global, List<BitSet> local) {
        int capacity = Math.max(16, global.size() * Math.max(1, local.size()));
        List<BitSet> candidates = new ArrayList<>(capacity);
        
        for (BitSet g : global) {
            for (BitSet t : local) {
                BitSet union = (BitSet) g.clone();
                union.or(t);
                candidates.add(union);
            }
        }
        
        return candidates;
    }
    
    /**
     * Runs the aggregation process for a given number of parts.
     * 
     * @param length Number of parts to process
     * @return List of global UCCs as BitSets
     * @throws IOException if file reading fails
     */
    public List<BitSet> runAggregation(int length) throws IOException {
        List<Path> partFiles = buildPartPaths(length);
        
        List<List<BitSet>> locals = new ArrayList<>();
        Integer attributeCount = null;
        
        for (Path file : partFiles) {
            UccFileReader.LocalUccResult result = readLocalUccFile(file);
            
            if (attributeCount == null) {
                attributeCount = result.attributeCount();
            } else if (attributeCount != result.attributeCount()) {
                throw new IllegalStateException(
                    "Attribute count mismatch: expected " + attributeCount + 
                    " but " + file + " has " + result.attributeCount()
                );
            }
            
            locals.add(result.uccs());
        }
        
        return aggregateGlobalUccs(locals, attributeCount);
    }
    
    /**
     * Prints the global UCCs to console.
     */
    public void printGlobalUccs(List<BitSet> global) {
        System.out.println("Global UCCs count = " + global.size());
        for (BitSet ucc : global) {
            System.out.println(UccMinimizer.toCommaList(ucc));
        }
    }
    
    /**
     * Simple static method to run aggregation with default settings.
     * Does everything in one call: reads files, aggregates, and prints results.
     * 
     * @param length Number of parts to process
     * @throws IOException if file reading fails
     */
    public static void run(int length) throws IOException {
        Path dir = Paths.get("serverReceived");
        String prefix = "serverpart_";
        String suffix = ".txt";

        UccAggregator agg = new UccAggregator(dir, prefix, suffix);
        List<BitSet> global = agg.runAggregation(length);
        agg.printGlobalUccs(global);
    }
//    public static void run(int length) throws IOException {
//        // 1) show where the JVM is running from
//        Path wd = Paths.get("").toAbsolutePath();
//        System.out.println("Working directory = " + wd);
//
//        // 2) resolve the folder you expect: <workingDir>/serverReceived
//        Path dir = wd.resolve("serverReceived").normalize().toAbsolutePath();
//        System.out.println("Input dir         = " + dir);
//
//        if (!java.nio.file.Files.isDirectory(dir)) {
//            throw new java.nio.file.NoSuchFileException(
//                    "serverReceived folder not found at: " + dir +
//                            "\nFix IntelliJ Run Configuration -> Working directory = C:/Users/jutt/IdeaProjects/ucc-bridge"
//            );
//        }
//
//        // 3) quick listing
//        System.out.println("Files in input dir:");
//        java.nio.file.Files.list(dir).forEach(p -> System.out.println("  - " + p.getFileName()));
//
//        String prefix = "serverpart_";
//        String suffix = ".txt";
//
//        UccAggregator agg = new UccAggregator(dir, prefix, suffix);
//
//        // 4) EXTRA DEBUG: verify each expected file exists
//        for (int i = 1; i <= length; i++) {
//            Path f = dir.resolve(prefix + i + suffix);
//            if (!java.nio.file.Files.isRegularFile(f)) {
//                throw new java.nio.file.NoSuchFileException("Missing file: " + f);
//            }
//        }
//
//        // 5) run + print
//        List<BitSet> global = agg.runAggregation(length);
//        agg.printGlobalUccs(global);
//    }

}
