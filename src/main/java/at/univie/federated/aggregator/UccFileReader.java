package at.univie.federated.aggregator;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.NoSuchFileException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;

/**
 * Handles reading and parsing UCC files.
 * Simple responsibility: read files and parse UCC data.
 */
public class UccFileReader {
    
    private final UccMinimizer minimizer;
    
    public UccFileReader() {
        this.minimizer = new UccMinimizer();
    }
    
    /**
     * Reads a local UCC file and returns parsed result.
     * 
     * Format:
     * - First line: number of attributes m (e.g., "9")
     * - Next lines: comma-separated attribute indexes (e.g., "1,5,6")
     */
    public LocalUccResult readLocalUccFile(Path file) throws IOException {
        if (!Files.isRegularFile(file)) {
            throw new NoSuchFileException("Missing part file: " + file);
        }
        
        int attributeCount;
        List<BitSet> uccs = new ArrayList<>();
        
        try (BufferedReader br = Files.newBufferedReader(file, StandardCharsets.UTF_8)) {
            String firstLine = br.readLine();
            if (firstLine == null) {
                throw new IOException("Empty file: " + file);
            }
            
            attributeCount = Integer.parseInt(firstLine.trim());
            
            String line;
            while ((line = br.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty()) {
                    continue;
                }
                
                BitSet ucc = parseUccLine(line, attributeCount, file);
                if (!ucc.isEmpty()) {
                    uccs.add(ucc);
                }
            }
        }
        
        // Defensively minimalize (local files should already be minimal)
        List<BitSet> minimal = minimizer.minimalize(uccs);
        return new LocalUccResult(attributeCount, minimal);
    }
    
    private BitSet parseUccLine(String line, int attributeCount, Path file) throws IOException {
        BitSet bs = new BitSet(attributeCount);
        
        for (String token : line.split(",")) {
            String trimmed = token.trim();
            if (trimmed.isEmpty()) {
                continue;
            }
            
            int index = Integer.parseInt(trimmed);
            if (index < 0 || index >= attributeCount) {
                throw new IOException("Attribute index out of range in " + file + ": " + index);
            }
            bs.set(index);
        }
        
        return bs;
    }
    
    /**
     * Data transfer object for local UCC result.
     * Java 11 compatible (replaces record which requires Java 14+).
     */
    public static class LocalUccResult {
        private final int attributeCount;
        private final List<BitSet> uccs;
        
        public LocalUccResult(int attributeCount, List<BitSet> uccs) {
            this.attributeCount = attributeCount;
            this.uccs = uccs;
        }
        
        public int attributeCount() {
            return attributeCount;
        }
        
        public List<BitSet> uccs() {
            return uccs;
        }
    }
}

