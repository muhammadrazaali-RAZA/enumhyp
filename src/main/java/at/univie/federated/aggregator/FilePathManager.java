package at.univie.federated.aggregator;

import java.io.IOException;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Objects;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Handles file path operations for UCC aggregator.
 * Simple responsibility: build and discover part file paths.
 */
public class FilePathManager {
    
    private final Path inputDir;
    private final String filePrefix;
    private final String fileSuffix;
    private final Pattern pattern;
    
    public FilePathManager(Path inputDir, String filePrefix, String fileSuffix) {
        this.inputDir = Objects.requireNonNull(inputDir);
        this.filePrefix = Objects.requireNonNull(filePrefix);
        this.fileSuffix = Objects.requireNonNull(fileSuffix);
        this.pattern = Pattern.compile(Pattern.quote(filePrefix) + "(\\d+)" + Pattern.quote(fileSuffix));
    }
    
    /**
     * Builds paths for known number of parts: serverpart_1.txt ... serverpart_n.txt
     */
    public List<Path> buildPartPaths(int n) {
        if (n <= 0) {
            throw new IllegalArgumentException("n must be > 0");
        }
        
        List<Path> paths = new ArrayList<>(n);
        for (int i = 1; i <= n; i++) {
            paths.add(inputDir.resolve(filePrefix + i + fileSuffix));
        }
        return paths;
    }
    
    /**
     * Auto-discovers files matching prefix+NUMBER+suffix and sorts by NUMBER.
     */
    public List<Path> discoverPartPaths() throws IOException {
        List<Path> paths = new ArrayList<>();
        
        try (DirectoryStream<Path> ds = Files.newDirectoryStream(inputDir)) {
            for (Path file : ds) {
                String name = file.getFileName().toString();
                Matcher m = pattern.matcher(name);
                if (m.matches()) {
                    paths.add(file);
                }
            }
        }
        
        paths.sort(Comparator.comparingInt(path -> extractIndex(path.getFileName().toString())));
        return paths;
    }
    
    private int extractIndex(String filename) {
        Matcher m = pattern.matcher(filename);
        if (!m.matches()) {
            return Integer.MAX_VALUE;
        }
        return Integer.parseInt(m.group(1));
    }
}

