package at.univie.federated.compared;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Compares two UCC result files and finds common data lines.
 * Simple responsibility: find and write common lines to output file.
 */
public class FileComparator {
    
    private final Path inputDir;
    private final String filePrefix;
    private final String fileSuffix;
    
    public FileComparator(Path inputDir, String filePrefix, String fileSuffix) {
        this.inputDir = inputDir;
        this.filePrefix = filePrefix;
        this.fileSuffix = fileSuffix;
    }
    
    /**
     * Compares multiple part files and finds lines common to ALL files.
     * 
     * @param length Number of files to compare (files from 1 to length)
     * @param outputFile Output file name (e.g., "common.txt")
     * @return Number of common lines found across all files
     * @throws IOException if file operations fail
     */
    public int compareAllAndWriteCommon(int length, String outputFile) throws IOException {
        if (length <= 0) {
            throw new IllegalArgumentException("Length must be > 0");
        }
        
        Path output = inputDir.resolve(outputFile);
        
        // Read first file to initialize common set
        Path firstFile = inputDir.resolve(filePrefix + 1 + fileSuffix);
        Set<String> commonLines = new HashSet<>(readDataLines(firstFile));
        
        // Intersect with remaining files
        for (int i = 2; i <= length; i++) {
            Path file = inputDir.resolve(filePrefix + i + fileSuffix);
            List<String> lines = readDataLines(file);
            Set<String> currentSet = new HashSet<>(lines);
            
            // Keep only lines that exist in both sets
            commonLines.retainAll(currentSet);
        }
        
        // Convert to sorted list for consistent output
        List<String> commonList = new java.util.ArrayList<>(commonLines);
        java.util.Collections.sort(commonList);
        
        // Write common lines to output file
        Files.write(output, commonList, StandardCharsets.UTF_8);
        
        return commonList.size();
    }
    
    
    /**
     * Reads data lines from a file, skipping the first line (attribute count).
     */
    private List<String> readDataLines(Path file) throws IOException {
        List<String> allLines = Files.readAllLines(file, StandardCharsets.UTF_8);
        
        // Skip first line (attribute count) and return data lines
        if (allLines.isEmpty()) {
            return new java.util.ArrayList<>();
        }
        
        return allLines.subList(1, allLines.size());
    }
    
    /**
     * Simple static method to compare all files with default settings.
     * Finds lines common across ALL files from 1 to length.
     * 
     * @param length Number of files to compare (files from 1 to length)
     * @throws IOException if file operations fail
     */
    public static void compare(int length) throws IOException {
        Path dir = Paths.get("serverReceived");
        String prefix = "serverpart_";
        String suffix = ".txt";
        
        FileComparator comparator = new FileComparator(dir, prefix, suffix);
        int commonCount = comparator.compareAllAndWriteCommon(length, "common.txt");
        
        System.out.println("Found " + commonCount + " common lines across all " + length + " files");
        System.out.println("Common lines written to: serverReceived/common.txt");
    }
}

